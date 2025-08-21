use tfhe::core_crypto::prelude::{polynomial_algorithms::polynomial_wrapping_add_assign, *};
use aligned_vec::ABox;
use concrete_fft::c64;
use crate::polys::*;

pub struct GaloisKeys {
    ggsw_fft: FourierGgswCiphertext<ABox<[c64]>>,
}


impl GaloisKeys {

    pub fn new_galois_keys(
        sk: &GlweSecretKeyOwned<u64>, 
        aut_idx: usize, 
        decomp_log: DecompositionBaseLog, 
        decomp_count: DecompositionLevelCount,
        distribution: Gaussian<f64>,
    ) -> Self {
        let mut seeder=  new_seeder();
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            seeder.seed(),
            seeder.as_mut(),
        );

        let aut_sk = sk
            .as_polynomial_list()
            .iter()
            .map(|poly| evaluate_automorphism(&poly, aut_idx))
            .collect::<Vec<_>>();

        let rest_base = 64 - decomp_log.0 * decomp_count.0;
        let mut ggsw = GgswCiphertext::new(
            0u64,
            sk.glwe_dimension().to_glwe_size(),
            sk.polynomial_size(),
            decomp_log,
            decomp_count,
            CiphertextModulus::new_native(),
        );

        encrypt_constant_ggsw_ciphertext(&sk, &mut ggsw, Plaintext(0), distribution, &mut encryption_generator);
        // correct noises
        ggsw.as_mut().iter_mut().for_each(|x| *x = (*x).wrapping_mul(1 << 16));

        let mut clearone = vec![0u64; sk.polynomial_size().0];
        clearone[0] = 1;
        for (level_idx, mut level_matrix) in ggsw.iter_mut().enumerate() {
            // idx from top to down
            let idx = decomp_count.0 - 1 - level_idx;

            // beta_j * m
            let mut clearpoly_beta_container = clearone.clone();
            clearpoly_beta_container[0] = clearpoly_beta_container[0].wrapping_mul(1 << (rest_base + decomp_log.0 * idx));
            let clearpoly_beta = Polynomial::from_container(clearpoly_beta_container);

            for (row_idx, mut row_as_glwe) in level_matrix.as_mut_glwe_list().iter_mut().enumerate()
            {
                let body = if row_idx + 1 < sk.glwe_dimension().to_glwe_size().0 {
                    let container = aut_sk[row_idx]
                        .as_ref()
                        .iter()
                        .map(|&v| {
                            v.wrapping_mul(1 << (rest_base + decomp_log.0 * idx)).wrapping_neg()
                        })
                        .collect::<Vec<_>>();

                    Polynomial::from_container(container)
                } else {
                    clearpoly_beta.clone()
                };

                // add it to the zero ciphertext `row_as_glwe`
                polynomial_wrapping_add_assign(&mut row_as_glwe.get_mut_body().as_mut_polynomial(), &body);
            }
        }

        let mut ggsw_fft = FourierGgswCiphertext::new(
            ggsw.glwe_size(), 
            ggsw.polynomial_size(), 
            ggsw.decomposition_base_log(), 
            ggsw.decomposition_level_count()
        );
        convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut ggsw_fft);

        Self { ggsw_fft }
    }

    pub fn keyswitch(&self, ct: &GlweCiphertextOwned<u64>) -> GlweCiphertextOwned<u64> {
        let mut ct_out = GlweCiphertext::new(
            0u64, 
            ct.glwe_size(), 
            ct.polynomial_size(), 
            CiphertextModulus::new_native()
        );

        add_external_product_assign(&mut ct_out, &self.ggsw_fft, &ct);
        ct_out
    }
}


#[test]
pub fn test_keyswitch() {
    let glwe_param = crate::params::DEFAULT_GLWE_PARAMTER;
    
    let mut seeder = new_seeder();
        
    let mut secret_generator =
        SecretRandomGenerator::<ActivatedRandomGenerator>::new(seeder.seed());

    let sk = GlweSecretKey::generate_new_binary(
        GlweDimension(glwe_param.glwe_size.0 - 1),
        glwe_param.polynomial_size,
        &mut secret_generator,
    );
    println!("sk len: {:?}", sk.as_polynomial_list().into_container().len());

    let mut seeder=  new_seeder();
    let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
        seeder.seed(),
        seeder.as_mut(),
    );

    let distribution=  Gaussian::from_dispersion_parameter(StandardDev(glwe_param.std_dev), 0.0);

    let cleartext = (0..glwe_param.polynomial_size.0).map(|i| {
        if i == 3 || i == 4 {
            glwe_param.delta
        } else {
            0
        }
    }).collect::<Vec<_>>();
    let pt = PlaintextList::from_container(cleartext);
    

    println!("gen ct");
    let mut ct = GlweCiphertext::new(0u64, glwe_param.glwe_size, glwe_param.polynomial_size, CiphertextModulus::new_native());
    encrypt_glwe_ciphertext(&sk, &mut ct, &pt, distribution, &mut encryption_generator);

    println!("gen galois keys");
    let aut_idx = glwe_param.polynomial_size.0 + 1;
    println!("aut idx: {}", aut_idx);
    let galois_keys = GaloisKeys::new_galois_keys(
        &sk, 
        aut_idx,
        glwe_param.decomposition_base_log, 
        glwe_param.decomposition_level_count, 
        distribution
    );

    println!("eval auto");
    let auto_ct = evaluate_automorphism_on_glwe_ciphertext(&ct, aut_idx);

    println!("ks");
    let switched_ct = galois_keys.keyswitch(&auto_ct);

    let mut out_pt = pt.clone();

    let auto_pt_poly = evaluate_automorphism(&pt.as_polynomial(), aut_idx);
    decrypt_glwe_ciphertext(&sk, &switched_ct, &mut out_pt);
    out_pt.as_mut().iter_mut().for_each(|v| *v /= glwe_param.delta / 16);
    println!("original pt: {:?}\n auto pt: {:?}\n dec pt: {:?}", pt, auto_pt_poly, out_pt);
}
