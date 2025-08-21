use std::collections::{HashMap, VecDeque};
use tfhe::core_crypto::prelude::polynomial_algorithms::*;
use tfhe::core_crypto::prelude::*;
use tfhe::core_crypto::algorithms::slice_algorithms::slice_wrapping_add_assign;
use crate::{params::*, keys::*};


pub struct Client {
    glwe_param: GlweParameter<u64>,
    pir_param: PirParam,

    sk: GlweSecretKeyOwned<u64>,
    distribution: Gaussian<f64>,
    threshold_size: usize,

    // PIR context
    row_depth: Vec<usize>,
    col_depth: Vec<usize>,

    // GGSW Database
    pub ppclient: ppverif_fhe::Client,
}

impl Client {
    pub fn new(glwe_param: GlweParameter<u64>, pir_param: PirParam) -> Self {
        assert!(glwe_param.polynomial_size.0 >= pir_param.row_num + pir_param.col_depth * glwe_param.decomposition_level_count.0);

        let mut seeder = new_seeder();
        
        let mut secret_generator =
            SecretRandomGenerator::<ActivatedRandomGenerator>::new(seeder.seed());

        let sk = GlweSecretKey::generate_new_binary(
            GlweDimension(glwe_param.glwe_size.0 - 1),
            glwe_param.polynomial_size,
            &mut secret_generator,
        );

        // calcualte row depth
        let mut row_depth = vec![0; pir_param.row_num];
        let mut queue = VecDeque::new();
        queue.push_back((pir_param.row_num, 0, 0));
        while let Some((rest, pos, depth)) = queue.pop_front() {
            if rest == 1 {
                row_depth[pos] = depth + 1;
                continue;
            }

            queue.push_back(((rest / 2) + (rest % 2), pos, depth + 1));
            queue.push_back((rest / 2, pos + (1 << depth), depth + 1));
        }

        // calculate column depth
        let mut col_depth = vec![0; pir_param.col_depth * glwe_param.decomposition_level_count.0];
        let mut queue = VecDeque::new();
        queue.push_back((col_depth.len(), 0, 0));
        while let Some((rest, pos, depth)) = queue.pop_front() {
            if rest == 1 {
                col_depth[pos] = depth + 1;
                continue;
            }

            queue.push_back(((rest / 2) + (rest % 2), pos, depth + 1));
            queue.push_back((rest / 2, pos + (1 << depth), depth + 1));
        }

        let ppclient = ppverif_fhe::Client::new(ppverif_fhe::DEFAULT_INNER_PRODUCT_PARAMETER, ppverif_fhe::DEFAULT_BLIND_ROTATION_PARAMETER);

        Self {
            glwe_param,
            pir_param,
            sk,
            threshold_size: glwe_param.polynomial_size.0 / 4,
            distribution: Gaussian::from_dispersion_parameter(StandardDev(glwe_param.std_dev), 0.0),
            row_depth,
            col_depth,
            ppclient,
        }
    }

    /// Build the galois keys corresponds to the possible automorphism
    /// In the current implementation, it is assumed that `row_num >= col_num`
    pub fn build_galois_keys(&self) -> HashMap<usize, GaloisKeys> {
        let mut map = HashMap::new();

        let mut row_queue = VecDeque::new();
        row_queue.push_back((self.pir_param.row_num, 1));

        let keys = GaloisKeys::new_galois_keys(
            &self.sk, 
            self.glwe_param.polynomial_size.0 + 1, 
            self.glwe_param.decomposition_base_log, 
            self.glwe_param.decomposition_level_count, 
            self.distribution
        );
        map.insert(self.glwe_param.polynomial_size.0 + 1, keys);

        while let Some((rest, depth)) = row_queue.pop_front() {
            if rest == 1 {
                continue;
            }

            let aut_idx = (self.glwe_param.polynomial_size.0 >> depth) + 1;
            let keys = GaloisKeys::new_galois_keys(
                &self.sk, 
                aut_idx, 
                self.glwe_param.decomposition_base_log, 
                self.glwe_param.decomposition_level_count, 
                self.distribution
            );

            map.insert(aut_idx, keys);

            row_queue.push_back(((rest / 2) + (rest % 2), depth + 1));
            row_queue.push_back((rest / 2, depth + 1));
        }

        map
    }

    /// The functional keyswtich keys are GGSW ciphertexts encrypting `-sk`
    pub fn build_functional_keyswitch_key(&self) -> Vec<GgswCiphertextOwned<u64>> {
        let mut seeder = new_seeder();
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            seeder.seed(),
            seeder.as_mut(),
        );

        let ciphertext_modulus = CiphertextModulus::new_native();

        // for each `s_i`, calcuate `GGSW(-s_i)`
        self.sk.as_polynomial_list()
            .iter()
            .map(|clearpoly| {
                let mut ggsw = GgswCiphertext::new(
                    0u64,
                    self.glwe_param.glwe_size,
                    self.glwe_param.polynomial_size,
                    self.glwe_param.decomposition_base_log,
                    self.glwe_param.decomposition_level_count,
                    ciphertext_modulus,
                );
        
                encrypt_constant_ggsw_ciphertext(&self.sk, &mut ggsw, Plaintext(0), self.distribution, &mut encryption_generator);
                ggsw.as_mut().iter_mut().for_each(|x| *x = (*x).wrapping_mul(self.glwe_param.factor));

                let rest_base = 64 - self.glwe_param.decomposition_base_log.0 * self.glwe_param.decomposition_level_count.0;
                for (level_idx, mut level_matrix) in ggsw.iter_mut().enumerate() {
                    // idx from top to down
                    let idx = self.glwe_param.decomposition_level_count.0 - 1 - level_idx;

                    // beta_j * m
                    let clearpoly_beta_container: Vec<_> = clearpoly
                        .iter()
                        .map(|vi| {
                            (*vi).wrapping_mul(1 << (rest_base + self.glwe_param.decomposition_base_log.0 * idx))
                        })
                        .collect();
                    let clearpoly_beta = Polynomial::from_container(clearpoly_beta_container);
        
                    for (row_idx, mut row_as_glwe) in level_matrix.as_mut_glwe_list().iter_mut().enumerate()
                    {
                        let mut body = Polynomial::from_container(vec![0; self.glwe_param.polynomial_size.0]);
        
                        if row_idx + 1 < self.glwe_param.glwe_size.0 {
                            // beta_j * s_i * m
                            polynomial_karatsuba_wrapping_mul(
                                &mut body,
                                &clearpoly_beta,
                                &self.sk.as_polynomial_list().get(row_idx),
                            );
        
                            // negate
                            body.as_mut().iter_mut().for_each(|v| *v = (*v).wrapping_neg());
                        } else {
                            // b
                            body.clone_from(&clearpoly_beta);
                        }
        
                        // subtract it to the zero ciphertext `row_as_glwe`, therefore the messages is `-sk`
                        polynomial_wrapping_sub_assign(&mut row_as_glwe.get_mut_body().as_mut_polynomial(), &body);
                    }
                }
        
        
                ggsw
            })
            .collect()
    }

    /// Construct the query GLWE ciphertext. The query idx should less than `row_num * col_num`.
    /// The query index is in column-first order.
    pub fn construct_query(&self, idx: usize) -> GlweCiphertextOwned<u64> {
        let row_idx = idx % self.pir_param.row_num;
        let col_idx = idx / self.pir_param.row_num;

        #[cfg(feature="verbose")]
        println!("row_idx: {}, col_idx: {}", row_idx, col_idx);
        let mut cleartext = vec![0; self.glwe_param.polynomial_size.0];
        
        cleartext[row_idx * 2] = self.glwe_param.delta >> self.row_depth[row_idx];
        let rest_base = 64 - self.glwe_param.decomposition_base_log.0 * self.glwe_param.decomposition_level_count.0;

        for depth in 0..self.pir_param.col_depth {
            let col_dig = (col_idx >> (self.pir_param.col_depth - 1 - depth)) & 0x1;
            #[cfg(feature="verbose")]
            println!("depth: {}, col_dig: {}", depth, col_dig);
            if col_dig == 1 {
                for beta_i in 0..self.glwe_param.decomposition_level_count.0 {
                    // mu * beta_i
                    let col_idx = depth * self.glwe_param.decomposition_level_count.0 + beta_i;
                    cleartext[
                        1 +  // column offset
                        2 * col_idx
                    ] = 1 << (rest_base + self.glwe_param.decomposition_base_log.0 * beta_i - self.col_depth[col_idx]);

                    #[cfg(feature="verbose")]
                    println!("col_idx: {}, mu: {}, depth: {}", col_idx, cleartext[
                        1 +  // column offset
                        2 * depth * self.glwe_param.decomposition_level_count.0 +   // group offset
                        2 * beta_i  // beta offset
                    ], self.col_depth[depth * self.glwe_param.decomposition_level_count.0 + beta_i]);
                }
            }
        }
        

        let mut ct = GlweCiphertext::new(0, self.glwe_param.glwe_size, self.glwe_param.polynomial_size, CiphertextModulus::new_native());
        let pt = PlaintextList::from_container(cleartext);
        
        let mut seeder=  new_seeder();
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            seeder.seed(),
            seeder.as_mut(),
        );

        encrypt_glwe_ciphertext(
            &self.sk, 
            &mut ct, 
            &PlaintextListOwned::new(0, pt.plaintext_count()), 
            self.distribution, 
            &mut encryption_generator,
        );

        ct.as_mut().iter_mut().for_each(|x| *x = (*x).wrapping_mul(self.glwe_param.factor));
        polynomial_wrapping_add_assign(&mut ct.get_mut_body().as_mut_polynomial(), &pt.as_polynomial());

        ct
    }

    /// GGSW(X^-a)
    pub fn construct_ggsw_monomial(&self, a: usize) -> GgswCiphertextOwned<u64> {
        let decomposition_base_log = DecompositionBaseLog(24);
        let decomposition_level_count = DecompositionLevelCount(1);

        let a = a % self.threshold_size;
        let mut ggsw = GgswCiphertext::new(
            0, 
            self.glwe_param.glwe_size, 
            self.glwe_param.polynomial_size, 
            decomposition_base_log, 
            decomposition_level_count, 
            CiphertextModulus::new_native()
        );

        let mut seeder=  new_seeder();
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            seeder.seed(),
            seeder.as_mut(),
        );
        
        encrypt_constant_ggsw_ciphertext(&self.sk, &mut ggsw, Plaintext(0), self.distribution, &mut encryption_generator);
        ggsw.as_mut().iter_mut().for_each(|x| *x = (*x).wrapping_mul(self.glwe_param.factor));

        let mut clearone = vec![0; self.glwe_param.polynomial_size.0];
        let poly_size = self.glwe_param.polynomial_size.0;
        clearone[(poly_size - a) % poly_size] = 1;
        // clearone[0] = 1;
        let clearpoly = Polynomial::from_container(clearone);

        let rest_idx = 64 - decomposition_base_log.0 * decomposition_level_count.0;
        for (level_idx, mut level_matrix) in ggsw.iter_mut().enumerate() {
            // idx from top to down
            let idx = decomposition_level_count.0 - 1 - level_idx;
            let mut clearpoly_beta = clearpoly.clone();
            // beta_j * m
            clearpoly_beta
                .as_mut()
                .iter_mut()
                .for_each(|vi| *vi <<= rest_idx + (decomposition_base_log.0 as usize * idx));

            for (row_idx, mut row_as_glwe) in level_matrix.as_mut_glwe_list().iter_mut().enumerate()
            {
                let mut body = Polynomial::from_container(vec![0; self.glwe_param.polynomial_size.0]);

                if row_idx + 1 < self.glwe_param.glwe_size.0 {
                    // beta_j * s_i * m
                    polynomial_wrapping_mul(
                        &mut body,
                        &clearpoly_beta,
                        &self.sk.as_polynomial_list().get(row_idx as usize),
                    );

                    // negate
                    body.as_mut().iter_mut().for_each(|v| *v = !(*v) + 1);
                } else {
                    // b
                    body.clone_from(&clearpoly_beta);
                }

                // add it to the zero ciphertext `row_as_glwe`
                slice_wrapping_add_assign(row_as_glwe.get_mut_body().as_mut(), body.as_ref());
            }
        }

        ggsw
    }

    pub fn construct_glwe_monomial(&self, a: usize) -> GlweCiphertextOwned<u64> {        
        let decomposition_base_log = 13;
        let a = a % self.threshold_size;
        let mut clearone = vec![0; self.glwe_param.polynomial_size.0];
        let poly_size = self.glwe_param.polynomial_size.0;
        clearone[(poly_size - a) % poly_size] = 1 << (64 - decomposition_base_log);
        let clearpoly = Polynomial::from_container(clearone);

        let mut glwe_ct = GlweCiphertext::new(
            0u64, 
            self.glwe_param.glwe_size, 
            self.glwe_param.polynomial_size, 
            CiphertextModulus::new_native()
        );
        
        let mut seeder=  new_seeder();
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
            seeder.seed(),
            seeder.as_mut(),
        );
        encrypt_glwe_ciphertext_assign(&self.sk, &mut glwe_ct, self.distribution, &mut encryption_generator);
        glwe_ct.as_mut().iter_mut().for_each(|x| *x = (*x).wrapping_mul(self.glwe_param.factor));

        polynomial_wrapping_add_assign(&mut glwe_ct.get_mut_body().as_mut_polynomial(), &clearpoly);

        glwe_ct
    }

    /// Decryption but under the delta to be `4 * glwe_params.delta` for debugging.
    pub fn decrypt(&self, ct: &GlweCiphertextOwned<u64>) -> Vec<u64> {
        let mut pt = PlaintextList::new(0, PlaintextCount(self.glwe_param.polynomial_size.0));
        decrypt_glwe_ciphertext(&self.sk, ct, &mut pt);

        pt.as_mut().iter_mut().for_each(|v| {
            let v_div = (*v as f64) / (self.glwe_param.delta as f64);
            *v = v_div.round() as u64 % self.glwe_param.plaintext_modulus;
        });
        pt.into_container()
    }

    pub fn decrypt_lwe(&self, ct: &LweCiphertextOwned<u64>) -> u64 {
        let pt = decrypt_lwe_ciphertext(&self.sk.as_lwe_secret_key(), ct);
        let v_div = pt.0 as f64 / ((self.glwe_param.delta >> 32) as f64);
        v_div.round() as u64 % self.glwe_param.plaintext_modulus
    }
}


#[test]
fn construct_client() {
    let client = Client::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER);
    println!("col depth: {:?}", client.col_depth);
    // assert_eq!(client.row_depth[0], 6);
}

#[test]
fn test_ggsw() {
    let mut glwe_param = DEFAULT_GLWE_PARAMTER;
    glwe_param.std_dev = 0.0;
    // glwe_param.decomposition_base_log.0 = 16;
    // glwe_param.decomposition_level_count.0 = 4;

    let mut pir_param = DEFAULT_PIR_PARAMETER;
    pir_param.row_num = 5;
    pir_param.col_depth = 1;
    // glwe_param.decomposition_level_count = DecompositionLevelCount(1);
    
    let mut seeder=  new_seeder();
    let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
        seeder.seed(),
        seeder.as_mut(),
    );

    let distribution=  Gaussian::from_dispersion_parameter(StandardDev(glwe_param.std_dev), 0.0);

    let client = Client::new(glwe_param, pir_param);
    let ksk = client.build_functional_keyswitch_key();
    let ksk_fft = ksk.iter().map(|kski| {
        let mut output_ggsw = FourierGgswCiphertext::new(
            glwe_param.glwe_size, 
            glwe_param.polynomial_size, 
            glwe_param.decomposition_base_log, 
            glwe_param.decomposition_level_count
        );
        convert_standard_ggsw_ciphertext_to_fourier(kski, &mut output_ggsw);
        output_ggsw
    }).collect::<Vec<_>>();

    let rest_base = 64 - glwe_param.decomposition_base_log.0 * glwe_param.decomposition_level_count.0;

    for (i, ggsw) in ksk_fft.iter().enumerate() {
        let mut cleartext = vec![0; glwe_param.polynomial_size.0];
        let idx = glwe_param.decomposition_level_count.0 - 1;
        cleartext[0] = 512 << (rest_base + idx * glwe_param.decomposition_base_log.0);
        
        let pt = PlaintextList::from_container(cleartext);
        let mut ct = GlweCiphertext::new(
            0, 
            glwe_param.glwe_size, 
            glwe_param.polynomial_size, 
            CiphertextModulus::new_native()
        );

        encrypt_glwe_ciphertext(&client.sk, &mut ct, &pt, distribution, &mut encryption_generator);

        let mut pt_out = pt.clone();
        polynomial_karatsuba_wrapping_mul(&mut pt_out.as_mut_polynomial(), &pt.as_polynomial(), &client.sk.as_polynomial_list().get(i));
        pt_out.as_mut().iter_mut().for_each(|v| *v >>= rest_base + idx * glwe_param.decomposition_base_log.0);
        println!("direct mul: {:?}", pt_out);


        let mut out = GlweCiphertext::new(
            0, 
            glwe_param.glwe_size, 
            glwe_param.polynomial_size, 
            CiphertextModulus::new_native()
        );
        let mut pt_dec = PlaintextList::new(0, PlaintextCount(glwe_param.polynomial_size.0));
        add_external_product_assign(&mut out, ggsw, &ct);
        decrypt_glwe_ciphertext(&client.sk, &out, &mut pt_dec);


        pt_dec.as_mut().iter_mut().for_each(|v| *v >>= rest_base + idx * glwe_param.decomposition_base_log.0);
        println!("dec: {:?},\n si: {:?}", pt_dec, client.sk.as_polynomial_list().get(i));
        break;
    }
}

#[test]
fn test_repack() {
    let mut glwe_param = DEFAULT_GLWE_PARAMTER;
    glwe_param.std_dev = 0.0;
    glwe_param.decomposition_base_log.0 = 10;
    glwe_param.decomposition_level_count.0 = 6;

    let mut pir_param = DEFAULT_PIR_PARAMETER;
    pir_param.row_num = 5;
    pir_param.col_depth = 1;
    // glwe_param.decomposition_level_count = DecompositionLevelCount(1);
    
    let mut seeder=  new_seeder();
    let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
        seeder.seed(),
        seeder.as_mut(),
    );

    let distribution=  Gaussian::from_dispersion_parameter(StandardDev(glwe_param.std_dev), 0.0);

    let client = Client::new(glwe_param, pir_param);
    let ksk = client.build_functional_keyswitch_key();
    let ksk_fft = ksk.iter().map(|kski| {
        let mut output_ggsw = FourierGgswCiphertext::new(
            glwe_param.glwe_size, 
            glwe_param.polynomial_size, 
            glwe_param.decomposition_base_log, 
            glwe_param.decomposition_level_count
        );
        convert_standard_ggsw_ciphertext_to_fourier(kski, &mut output_ggsw);
        output_ggsw
    }).collect::<Vec<_>>();


    // construct g = (1, beta, beta^2, ...)
    let rest_base = 64 - glwe_param.decomposition_base_log.0 * glwe_param.decomposition_level_count.0;
    let mut packed_ggsw = GgswCiphertext::new(
        0u64, 
        glwe_param.glwe_size, 
        glwe_param.polynomial_size, 
        glwe_param.decomposition_base_log, 
        glwe_param.decomposition_level_count, 
        CiphertextModulus::new_native()
    );
    
    for (level_idx, mut level_matrix) in packed_ggsw.iter_mut().enumerate() {
        // idx from top to down
        let idx = glwe_param.decomposition_level_count.0 - 1 - level_idx;
        let mut cleartext = vec![0; glwe_param.polynomial_size.0];
        
        cleartext[0] = 1 << (rest_base + idx * glwe_param.decomposition_base_log.0);

        // beta_j * m
        let pt = PlaintextList::from_container(cleartext);
        let mut ct = GlweCiphertext::new(
            0, 
            glwe_param.glwe_size, 
            glwe_param.polynomial_size, 
            CiphertextModulus::new_native()
        );
        encrypt_glwe_ciphertext(&client.sk, &mut ct, &pt, distribution, &mut encryption_generator);

        for (row_idx, mut row_as_glwe) in level_matrix.as_mut_glwe_list().iter_mut().enumerate()
        {
            if row_idx < ksk_fft.len() {
                add_external_product_assign(&mut row_as_glwe, &ksk_fft[row_idx], &ct);
            } else {
                glwe_ciphertext_add_assign(&mut row_as_glwe, &ct);
            }
        }
    }

    let dec_ggsw = decrypt_constant_ggsw_ciphertext(&client.sk, &packed_ggsw);
    println!("dec ggsw: {:?}", dec_ggsw);
    // should be ggsw encrypting a scalar

    let mut cleartext = vec![0; glwe_param.polynomial_size.0];
    let idx = glwe_param.decomposition_level_count.0 - 1;
    cleartext[0] = 512 << (rest_base + idx * glwe_param.decomposition_base_log.0);
    let pt = PlaintextList::from_container(cleartext);
    let mut ct = GlweCiphertext::new(
        0u64, 
        glwe_param.glwe_size, 
        glwe_param.polynomial_size, 
        CiphertextModulus::new_native()
    );

    polynomial_wrapping_add_assign(&mut ct.get_mut_body().as_mut_polynomial(), &pt.as_polynomial());
    encrypt_glwe_ciphertext_assign(&client.sk, &mut ct, distribution, &mut encryption_generator);
    let mut out = GlweCiphertext::new(
        0, 
        glwe_param.glwe_size, 
        glwe_param.polynomial_size, 
        CiphertextModulus::new_native()
    );
    let mut pt_dec = PlaintextList::new(0, PlaintextCount(glwe_param.polynomial_size.0));

    let mut packed_ggsw_fft = FourierGgswCiphertext::new(
        glwe_param.glwe_size, 
        glwe_param.polynomial_size, 
        glwe_param.decomposition_base_log, 
        glwe_param.decomposition_level_count
    );
    convert_standard_ggsw_ciphertext_to_fourier(&packed_ggsw, &mut packed_ggsw_fft);
    add_external_product_assign(&mut out, &packed_ggsw_fft, &ct);
    decrypt_glwe_ciphertext(&client.sk, &out, &mut pt_dec);
    pt_dec.as_mut().iter_mut().for_each(|v| *v >>= rest_base + idx * glwe_param.decomposition_base_log.0);
    println!("dec: {:?}", pt_dec);

}

#[test]
fn test_pir() {
    use rand::{thread_rng, RngCore};
    use crate::Server;

    let client = Client::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER);
    let galois_keys = client.build_galois_keys();
    let keyswitch_keys = client.build_functional_keyswitch_key();
    let server = Server::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER, galois_keys, keyswitch_keys);

    let mut rng = thread_rng();
    let now = std::time::Instant::now();
    // let query_idx = (1 << 15) + (1 << 14) + 1;
    let query_idx = rng.next_u64() as usize % (1 << 16);
    let query_ct = client.construct_query(query_idx);
    let elapsed = now.elapsed();
    println!("query time: {} millis", elapsed.as_micros() as f64 / 1000.0);

    
    let mut bs = (0..(1<<16)).map(|_| rng.next_u64() as usize % DEFAULT_GLWE_PARAMTER.polynomial_size.0).collect::<Vec<_>>();
    bs[query_idx] = 130;
    println!("starts counting");
    let now = std::time::Instant::now();
    let retrieved_ct = server.retrieve(&query_ct, &bs);
    let elapsed = now.elapsed();
    let pt = client.decrypt(&retrieved_ct);
    println!("pt: {:?}\n time: {} millis", pt, elapsed.as_micros() as f64 / 1000.0);
    println!("{}", pt.len());
}

#[test]
fn test_pir_with_ggsw() {
    use rand::{thread_rng, RngCore};
    use crate::Server;

    let client = Client::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER);
    let galois_keys = client.build_galois_keys();
    let keyswitch_keys = client.build_functional_keyswitch_key();
    let server = Server::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER, galois_keys, keyswitch_keys);

    let mut rng = thread_rng();
    // let query_idx = (1 << 15) + (1 << 14) + 1;
    let query_idx = rng.next_u64() as usize % (1 << 16);
    let query_ct = client.construct_query(query_idx);

    
    let mut bs = (0..(1<<16)).map(|_| rng.next_u64() as usize % DEFAULT_GLWE_PARAMTER.polynomial_size.0).collect::<Vec<_>>();
    bs[query_idx] = 130;
    let retrieved_ct = server.retrieve(&query_ct, &bs);


    let ggsw = client.construct_ggsw_monomial(2);
    let response = server.external_product(ggsw, retrieved_ct);

    let pt = client.decrypt(&response);
    println!("pt: {:?}", pt);
}

#[test]
fn test_pir_with_glwe() {
    use rand::{thread_rng, RngCore};
    use crate::Server;

    let client = Client::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER);
    let galois_keys = client.build_galois_keys();
    let keyswitch_keys = client.build_functional_keyswitch_key();
    let server = Server::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER, galois_keys, keyswitch_keys);

    let mut rng = thread_rng();
    // let query_idx = (1 << 15) + (1 << 14) + 1;
    let query_idx = rng.next_u64() as usize % (1 << 16);
    let query_ct = client.construct_query(query_idx);

    
    let mut bs = (0..(1<<16)).map(|_| rng.next_u64() as usize % DEFAULT_GLWE_PARAMTER.polynomial_size.0).collect::<Vec<_>>();
    bs[query_idx] = 130;
    let retrieved_ct = server.retrieve(&query_ct, &bs);


    let glwe = client.construct_glwe_monomial(2);
    let response = server.external_product_from_glwe(glwe, retrieved_ct);

    let pt = client.decrypt(&response);
    println!("pt: {:?}", pt);
}

#[test]
fn test_measure_client() {
    let mut client = Client::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER);
    
    let query_ct = client.construct_query(0);
    let query_len = query_ct.get_body().as_ref().len();  // body len
    let response_len = query_ct.get_mask().as_ref().len() + 1;  // LWE len

    let features = vec![0.0; 512];
    let template_ct = client.ppclient.encrypt_glwe(&features, 512.0);
    let template_len = template_ct.get_body().as_ref().len();  // for database size

    let ggsw = client.construct_ggsw_monomial(0);
    let ggsw_size = ggsw.as_glwe_list().iter().map(|glwe| glwe.get_body().as_ref().len()).sum::<usize>();
    let glwe_for_ggsw = client.construct_glwe_monomial(0);
    let glwe_for_ggsw_size = glwe_for_ggsw.get_body().as_ref().len();

    println!("query len: {}, template_len: {}, ggsw_size: {}, response_len: {}, glwe_for_ggsw: {}", query_len, template_len, ggsw_size, response_len, glwe_for_ggsw_size);
    println!(
        "query size: {} KB, template size: {} KB, ggsw size: {} KB, response size: {} KB, glwe_for_ggsw {} KB", 
        query_len as f32 * 8.0 / 1024.0, 
        template_len as f32 * 8.0 / 1024.0, 
        ggsw_size as f32 * 8.0 / 1024.0, 
        response_len as f32 * 4.0 / 1024.0,
        glwe_for_ggsw_size as f32 * 8.0 / 1024.0,
    );
}
