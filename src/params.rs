use tfhe::core_crypto::prelude::*;

#[derive(Debug, Copy, Clone)]
pub struct GlweParameter<Scalar> {
    pub glwe_size: GlweSize,
    pub polynomial_size: PolynomialSize,
    pub std_dev: f64,
    pub plaintext_modulus: Scalar,
    pub delta: Scalar,
    pub decomposition_base_log: DecompositionBaseLog,
    pub decomposition_level_count: DecompositionLevelCount,
    pub factor: Scalar,
}

pub const DEFAULT_GLWE_PARAMTER: GlweParameter<u64> = GlweParameter {
    glwe_size: GlweSize(5),
    polynomial_size: PolynomialSize(512),
    std_dev: 1.7347234759768072e-19,  // 3.2 / (2 ** 64)  // this will be multiplied by 2 ** 16 when encryption
    plaintext_modulus: 1 << 1,
    delta: 1 << (64 - 1),
    decomposition_base_log: DecompositionBaseLog(9),
    decomposition_level_count: DecompositionLevelCount(4),
    factor: 1 << 16,
};

#[derive(Debug, Copy, Clone)]
pub struct PirParam {
    pub row_num: usize,
    pub col_num: usize,
    pub col_depth: usize,
}

pub const DEFAULT_PIR_PARAMETER: PirParam = PirParam {
    row_num: 1 << 8,
    col_num: 1 << 8,
    col_depth: 8,
};