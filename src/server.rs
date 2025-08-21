use std::collections::{HashMap, VecDeque};
use tfhe::core_crypto::prelude::{polynomial_algorithms::polynomial_wrapping_monic_monomial_mul_assign, slice_algorithms::slice_wrapping_add_assign, *};
use crate::{keys::GaloisKeys, params::*, polys::*};

pub struct Server {
    glwe_param: GlweParameter<u64>,
    pir_param: PirParam,

    galois_keys: HashMap<usize, GaloisKeys>,
    keyswitch_keys: Vec<FourierGgswCiphertext<aligned_vec::ABox<[concrete_fft::c64]>>>,

    threshold_size: usize,
    lookup_table: PolynomialOwned<u64>,

    pub ppserver: ppverif_fhe::Server,
}


impl Server {
    pub fn new(
        glwe_param: GlweParameter<u64>, 
        pir_param: PirParam, 
        galois_keys: HashMap<usize, GaloisKeys>,
        keyswitch_keys: Vec<GgswCiphertextOwned<u64>>,
    ) -> Self {
        let keyswitch_keys = keyswitch_keys
            .into_iter()
            .map(|ksk| {
                let mut output_ggsw = FourierGgswCiphertext::new(
                    glwe_param.glwe_size, 
                    glwe_param.polynomial_size, 
                    glwe_param.decomposition_base_log, 
                    glwe_param.decomposition_level_count
                );
                convert_standard_ggsw_ciphertext_to_fourier(&ksk, &mut output_ggsw);
                output_ggsw
            })
            .collect();

        let threshold_size = glwe_param.polynomial_size.0 / 4;
        let ppserver = ppverif_fhe::Server::new(ppverif_fhe::DEFAULT_INNER_PRODUCT_PARAMETER, ppverif_fhe::DEFAULT_BLIND_ROTATION_PARAMETER, vec![]);
        Self {
            glwe_param,
            pir_param,
            galois_keys,
            keyswitch_keys,
            threshold_size,
            lookup_table: Self::build_lookup_table_from_threshold_size(threshold_size, glwe_param.polynomial_size),
            ppserver,
        }
    }

    /// Polynomial size `n` corresponds to the size of look-up table, 
    /// therefore the positive results are `{ idx in Z_n | idx > theta }`
    /// and the threshold size is the size of such set
    pub fn set_threshold_size(&mut self, threshold_size: usize) {
        assert!(threshold_size >= 1);
        self.threshold_size = threshold_size;
        self.lookup_table = Self::build_lookup_table_from_threshold_size(self.threshold_size, self.glwe_param.polynomial_size);
    }

    pub fn evaluate_automorphism_and_keyswitch(&self, ct: &GlweCiphertextOwned<u64>, aut_idx: usize) -> GlweCiphertextOwned<u64> {
        let aut_ct = evaluate_automorphism_on_glwe_ciphertext(ct, aut_idx);
        self.galois_keys.get(&aut_idx).unwrap().keyswitch(&aut_ct)
    }

    pub fn evaluate_automorphism_and_addassign(&self, ct: &mut GlweCiphertextOwned<u64>, aut_idx: usize) {
        let ct_aut = self.evaluate_automorphism_and_keyswitch(ct, aut_idx);
        glwe_ciphertext_add_assign(ct, &ct_aut);
    }

    pub fn act_query(&self, retrieved: &GlweCiphertextOwned<u64>, actor: &GgswCiphertextOwned<u64>) -> GlweCiphertextOwned<u64> {
        let mut out = GlweCiphertextOwned::new(
            0, 
            retrieved.glwe_size(), 
            retrieved.polynomial_size(), 
            CiphertextModulus::new_native()
        );

        let mut ggsw_fft = FourierGgswCiphertext::new(
            actor.glwe_size(), 
            actor.polynomial_size(), 
            actor.decomposition_base_log(), 
            actor.decomposition_level_count()
        );

        convert_standard_ggsw_ciphertext_to_fourier(&actor, &mut ggsw_fft);
        add_external_product_assign(&mut out, &ggsw_fft, &retrieved);

        out
    }

    pub fn retrieve(&self, query_ct: &GlweCiphertextOwned<u64>, bs: &[usize]) -> GlweCiphertextOwned<u64> {
        #[cfg(feature="verbose")]
        let now = std::time::Instant::now();
        let (glwe_rows, glwe_cols) = self.evaluate_trace(query_ct);
        #[cfg(feature="verbose")]
        println!("trace time: {} millis", now.elapsed().as_micros() as f64 / 1000.0);

        #[cfg(feature="verbose")]
        let now = std::time::Instant::now();
        let ggsw_cols = glwe_cols
            .chunks(self.glwe_param.decomposition_level_count.0)
            .map(|col_cts| {
                let ggsw = self.glwe_ciphertexts_to_ggsw_ciphertext(col_cts);
                let mut ggsw_fft = FourierGgswCiphertext::new(
                    self.glwe_param.glwe_size, 
                    self.glwe_param.polynomial_size, 
                    self.glwe_param.decomposition_base_log, 
                    self.glwe_param.decomposition_level_count
                );
                convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut ggsw_fft);
                ggsw_fft
            })
            .collect::<Vec<_>>();
        #[cfg(feature="verbose")]
        println!("pack + fft time: {} millis", now.elapsed().as_micros() as f64 / 1000.0);

        // column-first bs
        #[cfg(feature="verbose")]
        let now = std::time::Instant::now();
        let row_expansion = glwe_rows
            .into_iter()
            .map(|row_ct| self.expand_glwe_ciphertext_multiplying_lookup_table(&row_ct))
            .collect::<Vec<_>>();
        #[cfg(feature="verbose")]
        println!("expansion time: {} millis", now.elapsed().as_micros() as f64 / 1000.0);

        #[cfg(feature="verbose")]
        let now = std::time::Instant::now();
        let mut row_folded = bs
            .chunks(self.pir_param.row_num)
            .map(|bi| self.glwe_ciphertexts_innerprod_along_columns(&row_expansion, bi))
            .collect::<Vec<_>>();
        #[cfg(feature="verbose")]
        println!("row folding time: {} millis", now.elapsed().as_micros() as f64 / 1000.0);

        #[cfg(feature="verbose")]
        let now = std::time::Instant::now();
        let mut row_folded_slice = &mut row_folded[..];
        for ggsw_fft in ggsw_cols {
            let len = row_folded_slice.len();
            let (lefts, rights) = row_folded_slice.split_at_mut(len / 2);

            for (ct_left, ct_right) in lefts.iter_mut().zip(rights.iter_mut()) {
                cmux_assign(ct_left, ct_right, &ggsw_fft);
            }

            row_folded_slice = &mut row_folded[..len/2];
        }
        #[cfg(feature="verbose")]
        println!("col folding time: {} millis", now.elapsed().as_micros() as f64 / 1000.0);

        row_folded[0].clone()
    }

    fn build_lookup_table_from_threshold_size(threshold_size: usize, poly_size: PolynomialSize) -> PolynomialOwned<u64> {
        let mut container = vec![0; poly_size.0];
        container[(poly_size.0/2 - threshold_size)..poly_size.0/2]
            .iter_mut()
            .for_each(|v| *v = 1);

        Polynomial::from_container(container)
    }

    fn expand_glwe_ciphertext_multiplying_lookup_table(&self, glwe_ct: &GlweCiphertextOwned<u64>) -> Vec<Vec<u64>> {
        glwe_ct.as_polynomial_list()
            .iter()
            .map(|poly| {
                let n = poly.polynomial_size().0;

                let mut mul_lut = Polynomial::new(0, poly.polynomial_size());

                let mut conv_coe = vec![0; n + (n / 2) - 1];
                conv_coe[(n/2)-1..].copy_from_slice(poly.as_ref());
                for (tgt, &src) in conv_coe[..(n/2)-1].iter_mut().zip(poly.as_ref()[(n/2)+1..].iter()) {
                    *tgt = src.wrapping_neg();
                }

                mul_lut[0] = conv_coe[..self.threshold_size].iter().fold(0, |acc, &v| acc.wrapping_add(v));
                for i in 1..n {
                    mul_lut[i] = mul_lut[i - 1].wrapping_add(conv_coe[i + self.threshold_size - 1]);
                    mul_lut[i] = mul_lut[i].wrapping_sub(conv_coe[i - 1]);
                }
                // same to the following function but is with only O(n) complexity
                // polynomial_karatsuba_wrapping_mul(&mut mul_lut, &poly, &self.lookup_table);

                let mut conv_res = vec![0; n + n - 1];
                conv_res[..n].copy_from_slice(mul_lut.as_ref());
                for i in n..(n+n-1) {
                    conv_res[i] = conv_res[i - n].wrapping_neg();
                }

                conv_res
            })
            .collect()
    }

    fn glwe_ciphertexts_innerprod_along_columns(&self, expanded_cts: &Vec<Vec<Vec<u64>>>, bs: &[usize]) -> GlweCiphertextOwned<u64> {
        let mut ct = GlweCiphertext::new(
            0u64, 
            self.glwe_param.glwe_size, 
            self.glwe_param.polynomial_size, 
            CiphertextModulus::new_native()
        );

        let n = self.glwe_param.polynomial_size.0;
        for (expanded_ct, &b) in expanded_cts.iter().zip(bs.iter()) {
            let b = b % n; 
            for (expanded_poly, mut ct_i) in expanded_ct.iter().zip(ct.as_mut_polynomial_list().iter_mut()) {
                slice_wrapping_add_assign(&mut ct_i.as_mut(), &expanded_poly[b..b+n]);
            }
        }

        ct
    }

    fn glwe_ciphertext_wrapping_monic_monomial_mul_assign(&self, ct: &mut GlweCiphertextOwned<u64>, monomial_degree: MonomialDegree) {
        for mut ct_i in ct.as_mut_polynomial_list().iter_mut() {
            polynomial_wrapping_monic_monomial_mul_assign(&mut ct_i, monomial_degree);
        }
    }

    /// Given a packed ciphertext, extracte the coefficients to GLWE ciphertexts
    fn evaluate_trace(&self, ct: &GlweCiphertextOwned<u64>) -> (Vec<GlweCiphertextOwned<u64>>, Vec<GlweCiphertextOwned<u64>>) {
        let mut row_glwes = vec![None; self.pir_param.row_num];
        let mut col_glwes = vec![None; self.pir_param.col_depth * self.glwe_param.decomposition_level_count.0];
        // let mut col_glwes = vec![None; self.pir_param.col_depth];

        // split row and col cts
        let mut ct_row = ct.clone();
        self.evaluate_automorphism_and_addassign(&mut ct_row, self.glwe_param.polynomial_size.0 + 1);

        // =========== Logically, eval permutation then add ==============
        // let mut ct_col = ct.clone();
        // self.glwe_ciphertext_wrapping_monic_monomial_mul_assign(&mut ct_col, MonomialDegree(2 * self.glwe_param.polynomial_size.0 - 1));
        // self.evaluate_automorphism_and_addassign(&mut ct_col, self.glwe_param.polynomial_size.0 + 1);
        // =========== Can direct sub extracted then rotate ==============
        let ct_col_cont = ct.as_ref().iter().zip(ct_row.as_ref().iter()).map(|(&cti, &ct_rowi)| {
            cti.wrapping_mul(2).wrapping_sub(ct_rowi)
        }).collect::<Vec<_>>();
        let mut ct_col = GlweCiphertext::from_container(ct_col_cont, ct.polynomial_size(), ct.ciphertext_modulus());
        self.glwe_ciphertext_wrapping_monic_monomial_mul_assign(&mut ct_col, MonomialDegree(2 * self.glwe_param.polynomial_size.0 - 1));

        // handle rows
        let mut row_queue = VecDeque::new();
        row_queue.push_back((ct_row, self.pir_param.row_num, 0, 1));

        while let Some((mut ct_i, rest, pos, depth)) = row_queue.pop_front() {
            if rest == 1 {
                row_glwes[pos / 2] = Some(ct_i);
                continue;
            }
            
            // automoprhism without shifting
            let mut ct_i_no_shift = ct_i.clone();
            self.evaluate_automorphism_and_addassign(&mut ct_i_no_shift, (self.glwe_param.polynomial_size.0 >> depth) + 1);
            let even_rest = (rest / 2) + (rest % 2);
            
            
            // automorphism with shifting
            ct_i.as_mut().iter_mut().zip(ct_i_no_shift.as_ref().iter()).for_each(|(shift_i, no_shift_i)| {
                *shift_i = shift_i.wrapping_mul(2).wrapping_sub(*no_shift_i);
            });
            self.glwe_ciphertext_wrapping_monic_monomial_mul_assign(&mut ct_i, MonomialDegree(2 * self.glwe_param.polynomial_size.0 - (1 << depth)));
            // Also replace automorphism to "sub-then-rotate"
            // self.evaluate_automorphism_and_addassign(&mut ct_i, (self.glwe_param.polynomial_size.0 >> depth) + 1);

            row_queue.push_back((ct_i_no_shift, even_rest, pos, depth + 1));
            row_queue.push_back((ct_i, rest / 2, pos + (1 << depth), depth + 1));
        }

        // handle columns
        let mut col_queue = VecDeque::new();
        col_queue.push_back((ct_col, col_glwes.len(), 0, 1));

        while let Some((mut ct_i, rest, pos, depth)) = col_queue.pop_front() {
            // println!("eval trace: {}, {}, {}", rest, pos, depth);
            if rest == 1 {
                // println!("done");
                col_glwes[pos / 2] = Some(ct_i);
                continue;
            }
            

            // automoprhism without shifting
            let mut ct_i_no_shift = ct_i.clone();
            self.evaluate_automorphism_and_addassign(&mut ct_i_no_shift, (self.glwe_param.polynomial_size.0 >> depth) + 1);
            let even_rest = (rest / 2) + (rest % 2);

            
            // automorphism with shifting
            ct_i.as_mut().iter_mut().zip(ct_i_no_shift.as_ref().iter()).for_each(|(shift_i, no_shift_i)| {
                *shift_i = shift_i.wrapping_mul(2).wrapping_sub(*no_shift_i);
            });
            self.glwe_ciphertext_wrapping_monic_monomial_mul_assign(&mut ct_i, MonomialDegree(2 * self.glwe_param.polynomial_size.0 - (1 << depth)));
            // self.evaluate_automorphism_and_addassign(&mut ct_i, (self.glwe_param.polynomial_size.0 >> depth) + 1);

            col_queue.push_back((ct_i_no_shift, even_rest, pos, depth + 1));
            col_queue.push_back((ct_i, rest / 2, pos + (1 << depth), depth + 1));
        }

        let row_glwes = row_glwes
            .into_iter()
            .map(|some_ct| some_ct.unwrap())
            .collect();
        let col_glwes = col_glwes
            .into_iter()
            .map(|some_ct| some_ct.unwrap())
            .collect();

        (row_glwes, col_glwes)
    }

    /// Pack glwe ciphertexts to a single ggsw ciphertext.
    /// The glwe ciphertexts should be sorted to encrypt `m * g`
    fn glwe_ciphertexts_to_ggsw_ciphertext(&self, glwe_cts: &[GlweCiphertextOwned<u64>]) -> GgswCiphertextOwned<u64> {
        let mut ggsw = GgswCiphertext::new(
            0u64, 
            self.glwe_param.glwe_size, 
            self.glwe_param.polynomial_size, 
            self.glwe_param.decomposition_base_log, 
            self.glwe_param.decomposition_level_count, 
            CiphertextModulus::new_native()
        );

        for (mut level_matrix, glwe_ct) in ggsw.iter_mut().zip(glwe_cts.iter().rev()) {
            for (row_idx, mut row_as_glwe) in level_matrix.as_mut_glwe_list().iter_mut().enumerate()
            {
                if row_idx < self.keyswitch_keys.len() {
                    add_external_product_assign(&mut row_as_glwe, &self.keyswitch_keys[row_idx], &glwe_ct);
                } else {
                    glwe_ciphertext_add_assign(&mut row_as_glwe, &glwe_ct);
                }
            }
        }

        ggsw
    }

    pub fn external_product(&self, ggsw: GgswCiphertextOwned<u64>, glwe: GlweCiphertextOwned<u64>) -> GlweCiphertextOwned<u64> {
        let mut ggsw_fft = FourierGgswCiphertext::new(
            ggsw.glwe_size(), 
            ggsw.polynomial_size(), 
            ggsw.decomposition_base_log(), 
            ggsw.decomposition_level_count()
        );

        convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut ggsw_fft);
        let mut out = GlweCiphertext::new(
            0, 
            glwe.glwe_size(), 
            glwe.polynomial_size(), 
            glwe.ciphertext_modulus()
        );
        add_external_product_assign(&mut out, &ggsw_fft, &glwe);
        out
    }

    pub fn external_product_from_glwe(&self, glwe_to_ggsw: GlweCiphertextOwned<u64>, glwe: GlweCiphertextOwned<u64>) -> GlweCiphertextOwned<u64> {
        // glwe to ggsw
        let mut ggsw = GgswCiphertext::new(
            0u64, 
            self.glwe_param.glwe_size, 
            self.glwe_param.polynomial_size, 
            DecompositionBaseLog(13), 
            DecompositionLevelCount(1), 
            CiphertextModulus::new_native()
        );

        for (mut level_matrix, glwe_ct) in ggsw.iter_mut().zip([glwe_to_ggsw].into_iter()) {
            for (row_idx, mut row_as_glwe) in level_matrix.as_mut_glwe_list().iter_mut().enumerate()
            {
                if row_idx < self.keyswitch_keys.len() {
                    add_external_product_assign(&mut row_as_glwe, &self.keyswitch_keys[row_idx], &glwe_ct);
                } else {
                    glwe_ciphertext_add_assign(&mut row_as_glwe, &glwe_ct);
                }
            }
        }

        self.external_product(ggsw, glwe)
    }
}


#[test]
fn test_extraction() {
    let mut glwe_param = DEFAULT_GLWE_PARAMTER;
    glwe_param.std_dev = 0.0;
    let mut pir_param = DEFAULT_PIR_PARAMETER;
    pir_param.row_num = 5;
    pir_param.col_depth = 1;
    // glwe_param.decomposition_level_count = DecompositionLevelCount(1);
    
    let mut seeder = new_seeder();
        
    let mut secret_generator =
        SecretRandomGenerator::<ActivatedRandomGenerator>::new(seeder.seed());

    let sk = GlweSecretKey::generate_new_binary(
        GlweDimension(glwe_param.glwe_size.0 - 1),
        glwe_param.polynomial_size,
        &mut secret_generator,
    );

    let mut seeder=  new_seeder();
    let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
        seeder.seed(),
        seeder.as_mut(),
    );

    let distribution=  Gaussian::from_dispersion_parameter(StandardDev(glwe_param.std_dev), 0.0);

    let cleartext = (0..glwe_param.polynomial_size.0).map(|i| {
        if i <= 9 && i % 2 == 0 {
            i as u64 * glwe_param.delta / 16
        } else if i <= 17 {
            (i as u64 / 2) * glwe_param.delta / 16
        } else {
            0
        }
    }).collect::<Vec<_>>();
    let pt = PlaintextList::from_container(cleartext);
    
    let mut ct = GlweCiphertext::new(0u64, glwe_param.glwe_size, glwe_param.polynomial_size, CiphertextModulus::new_native());
    encrypt_glwe_ciphertext(&sk, &mut ct, &pt, distribution, &mut encryption_generator);

    let galois_keys = (0..=5).fold(HashMap::new(), |mut map, i| {
        let keys=  GaloisKeys::new_galois_keys(
            &sk, 
            (glwe_param.polynomial_size.0 >> i) + 1,
            glwe_param.decomposition_base_log, 
            glwe_param.decomposition_level_count, 
            distribution
        );

        map.insert((glwe_param.polynomial_size.0 >> i) + 1, keys);
        map
    });

    let server = Server::new(glwe_param, pir_param, galois_keys, vec![]);
    let (row_cts, _col_cts) = server.evaluate_trace(&ct);
    for (i, ct_i) in row_cts.iter().enumerate() {
        let mut pt = pt.clone();
        decrypt_glwe_ciphertext(&sk, &ct_i, &mut pt);
        pt.as_mut().iter_mut().for_each(|v| *v /= glwe_param.delta / 16);
        println!("i: {}, dec_pt: {:?}", i, pt);
    }
}

#[test]
fn test_poly() {
    use rand::{thread_rng, RngCore};
    use tfhe::core_crypto::prelude::polynomial_algorithms::*;
    let server = Server::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER, HashMap::new(), vec![]);

    let mut glwe_ct = GlweCiphertextOwned::new(0u64, DEFAULT_GLWE_PARAMTER.glwe_size, DEFAULT_GLWE_PARAMTER.polynomial_size, CiphertextModulus::new_native());
    let mut rng = thread_rng();
    glwe_ct.as_mut().iter_mut().for_each(|v| *v = rng.next_u64());

    println!("expansion");
    let expanded = server.expand_glwe_ciphertext_multiplying_lookup_table(&glwe_ct);
    let n = DEFAULT_GLWE_PARAMTER.polynomial_size.0;
    for (expansion, ct_i) in expanded.iter().zip(glwe_ct.as_polynomial_list().iter()) {
        for j in 1..n {
            let lut_mul = server.lookup_table.clone();
            // println!("lut: {:?}", lut_mul);
            let mut lut_out = Polynomial::new(0, lut_mul.polynomial_size());

            if j > 0 {
                polynomial_wrapping_monic_monomial_mul(&mut lut_out, &server.lookup_table, MonomialDegree(n + n - j));
            }
            // println!("muled: {:?}", lut_out);

            let mut lut_mul_ct = lut_out.clone();
            polynomial_karatsuba_wrapping_mul(&mut lut_mul_ct, &lut_out, &ct_i);

            let lut_cons = Polynomial::from_container(&expansion[j..j+n]);
            polynomial_wrapping_sub_assign(&mut lut_mul_ct, &lut_cons);

            // println!("diff: {:?}", lut_mul_ct);
            assert_eq!(lut_mul_ct.into_container(), vec![0; n]);
            // break;
        }
        break;
    }
}
