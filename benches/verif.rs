use criterion::{criterion_group, criterion_main, Criterion};
use anonyverif::*;

fn anonymous_verification_benchmark(c: &mut Criterion) {
    let mut client = Client::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER);
    let features = vec![0.0; 512];
    c.bench_function("Encrypt Template", |b| b.iter(|| client.ppclient.encrypt_glwe(&features, 512.0)));


    for (row_num, col_num, col_depth) in [
        (1 << 5, 1 << 5, 5), // 2 ** 10
        (1 << 6, 1 << 5, 5), // 2 ** 11
        (1 << 5, 1 << 6, 6),
        (1 << 6, 1 << 6, 6), // 2 ** 12
        (1 << 7, 1 << 5, 5), 
        (1 << 8, 1 << 4, 4),
        (1 << 7, 1 << 6, 6), // 2 ** 13
        (1 << 8, 1 << 5, 5),
        (1 << 7, 1 << 7, 7), // 2 ** 14
        (1 << 8, 1 << 6, 6),
        (1 << 8, 1 << 7, 7), // 2 ** 15
        (1 << 8, 1 << 8, 8), // 2 ** 16
    ] {
        let client = Client::new(DEFAULT_GLWE_PARAMTER, PirParam { row_num, col_num, col_depth });
        c.bench_function(&format!("Build Query {}/row={}/col={}", row_num * col_num, row_num, col_num), |b| 
            b.iter(|| client.construct_query(0))
        );

        let galois_keys = client.build_galois_keys();
        let keyswitch_keys = client.build_functional_keyswitch_key();
        let server = Server::new(DEFAULT_GLWE_PARAMTER, PirParam { row_num, col_num, col_depth }, galois_keys, keyswitch_keys);

        let bs = vec![0; row_num * col_num];
        let query_ct = client.construct_query(0);
        c.bench_function(&format!("PIR Response {}/row={}/col={}", row_num * col_num, row_num, col_num), |b| {
            b.iter(|| server.retrieve(&query_ct, &bs));
        });
    }
}

criterion_group!(benches, anonymous_verification_benchmark);
criterion_main!(benches);