# Anonymous Verification

This repository provides code example and benchmark test on the anonymous verification scheme.
The dependency repository [ppverif_fhe](https://github.com/Fainabi/ppverif_fhe) is the underline similarity implementation. 

## Dependencies

Rust release version >= 1.84 to compile TFHE-rs.


## Example

We provide a demo for this repository:
```sh
$ cargo run --release --example demo
```


To run the benchmark, type:
```sh
$ cargo bench --bench matmul
```

We also provide `block` and `verif` bench for the building blocks and PIR benchmark.
