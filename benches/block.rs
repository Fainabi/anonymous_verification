use criterion::{criterion_group, criterion_main, Criterion};
use anonyverif::*;
use tfhe::core_crypto::prelude::*;

fn building_block_benchmark(c: &mut Criterion) {
    let client = Client::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER);
    let galois_keys = client.build_galois_keys();
    let autidx = galois_keys.keys().next().copied().unwrap();
    let keyswitch_keys = client.build_functional_keyswitch_key();
    let server = Server::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER, galois_keys, keyswitch_keys);

    let query_ct = client.construct_query(0);

    c.bench_function("Automorphism", |b| b.iter(|| {
        server.evaluate_automorphism_and_keyswitch(&query_ct, autidx);
    }));
    
    let mut ct3 = query_ct.clone();
    c.bench_function("HomAdd", |b| b.iter(|| {
        glwe_ciphertext_add_assign(&mut ct3, &query_ct);
    }));
    
    let mut ct = query_ct.clone();
    let mut seeder = new_seeder();
    let mut secret_generator =
        SecretRandomGenerator::<ActivatedRandomGenerator>::new(seeder.seed());
    let sk = GlweSecretKey::generate_new_binary(
        GlweDimension(DEFAULT_GLWE_PARAMTER.glwe_size.0 - 1),
        DEFAULT_GLWE_PARAMTER.polynomial_size,
        &mut secret_generator,
    );
    let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(
        seeder.seed(),
        seeder.as_mut(),
    );

    let mut ggsw = GgswCiphertext::new(0u64, DEFAULT_GLWE_PARAMTER.glwe_size, DEFAULT_GLWE_PARAMTER.polynomial_size, DEFAULT_GLWE_PARAMTER.decomposition_base_log, DEFAULT_GLWE_PARAMTER.decomposition_level_count, CiphertextModulus::new_native());
    encrypt_constant_ggsw_ciphertext(&sk, &mut ggsw, Plaintext(0), Gaussian::from_dispersion_parameter(StandardDev(3.2), 0.0), &mut encryption_generator);
    let mut ggsw_fft = FourierGgswCiphertext::new(ggsw.glwe_size(), ggsw.polynomial_size(), ggsw.decomposition_base_log(), ggsw.decomposition_level_count());
    
    c.bench_function("ExternalProduct", |b| b.iter(|| {
        convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut ggsw_fft);
        add_external_product_assign(&mut ct, &ggsw_fft, &query_ct);
    }));

    c.bench_function("GGSW Encryption", |b| b.iter(|| {
        client.construct_ggsw_monomial(0);
    }));
}

criterion_group!(benches, building_block_benchmark);
criterion_main!(benches);