use tfhe::core_crypto::prelude::*;

pub(crate) fn evaluate_automorphism<Scalar, Cont>(
    poly: &Polynomial<Cont>,
    aut_idx: usize,
) -> PolynomialOwned<Scalar>
where
    Scalar: UnsignedInteger,
    Cont: Container<Element = Scalar>,
{
    let poly_size = poly.polynomial_size().0;
    let mut out_cont = vec![Scalar::ZERO; poly_size];
    
    for (i, &v) in poly.iter().enumerate() {
        let new_idx = (i * aut_idx) % (2 * poly_size);
        if new_idx < poly_size {
            out_cont[new_idx] = out_cont[new_idx].wrapping_add(v);
        } else {
            out_cont[new_idx - poly_size] = out_cont[new_idx - poly_size].wrapping_sub(v);
        }
    }

    Polynomial::from_container(out_cont)
}

pub(crate) fn evaluate_automorphism_on_glwe_ciphertext<Scalar, Cont>(
    ct: &GlweCiphertext<Cont>,
    aut_idx: usize
) -> GlweCiphertextOwned<Scalar>
where
    Scalar: UnsignedInteger,
    Cont: Container<Element = Scalar>,
{
    let container = ct
        .as_polynomial_list()
        .iter()
        .flat_map(|poly| evaluate_automorphism(&poly, aut_idx).into_container().into_iter())
        .collect::<Vec<_>>();

    GlweCiphertext::from_container(container, ct.polynomial_size(), CiphertextModulus::new_native())
}
