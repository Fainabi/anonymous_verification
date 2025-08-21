use criterion::{criterion_group, criterion_main, Criterion};
use anonyverif::*;
use rand::*;

fn matmul_benchmark(c: &mut Criterion) {
    let mut client = Client::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER);
    client.ppclient.enroll_ggsw_masks(0);

    let mat = client.ppclient.collect_database_into_matrix().unwrap();
    let ncol = mat.column_iter().count();
    // let ncol = 512 * 5;
    let mut rng = thread_rng();
    let col_vec = nalgebra::DVector::from_iterator(ncol, (0..ncol).map(|_| rng.next_u64()));
    for nrow_log in 10..16 {
        let nrow = 1usize << nrow_log;
        let mat = nalgebra::DMatrix::from_iterator(nrow, ncol, (0..nrow*ncol).map(|_| rng.next_u64()));

        c.bench_function(&format!("nrow={}", nrow), |b| {
            b.iter(|| {
                let _ = &mat * &col_vec;
            });
        });
    }
}

criterion_group!(benches, matmul_benchmark);
criterion_main!(benches);