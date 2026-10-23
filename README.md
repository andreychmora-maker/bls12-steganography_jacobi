# Optimizing Elligator Squared for BLS Curves via Jacobi Quartics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23537349.svg)](https://doi.org/10.5281/zenodo.23537349)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the Magma Computer Algebra System proof-of-concept implementation for the steganographic obfuscation of BLS12-family public keys. 

Specifically, this repository demonstrates a highly optimized Elligator Squared pipeline for the **BLS12-479+** and **BLS12-539+** curves. By leveraging their algebraically engineered 2-torsion subgroups, the pipeline constructs a mathematically exact birational equivalence to the **Extended Jacobi Quartic** space, effectively shattering the traditional computational penalty of public key obfuscation.

## Mathematical & Architectural Highlights

1. **Inversion-Free Arithmetic**: By carefully fusing the Billet-Joye birational mappings with Hişil's homogeneous projective coordinates $(X:Y:Z)$, we completely eliminate intermediate field inversions. The pipeline strictly defers the solitary field inversion to the final affine projection at the very end of the deobfuscation cycle.
2. **Jacobi Quartic Sandbox**: Heavy cryptographic operations are executed entirely within the accelerated Jacobi space, utilizing an optimized dedicated doubling formula ($3\mathbf{M}+4\mathbf{S}+1\mathbf{C}+3\mathbf{add}$).
3. **Dual Operational Profiles**: 
   * **High-Throughput Mode (Asymmetric $w\text{NAF}$)**: Employs a wide window ($w=5$) for full-scalar blinding and a minimal window ($w=2$) for cofactor clearing. Precomputation tables are normalized ($Z=1$) using Montgomery Batch Inversion, enabling the exclusive use of dedicated addition formulas ($9\mathbf{M}+5\mathbf{S}+1\mathbf{C}+6\mathbf{add}$).
   * **SPA-Resistant Constant-Time Mode (Montgomery Ladder)**: Replaces $w\text{NAF}$ with a strictly constant-time Montgomery ladder over the Jacobi Quartic to thwart Simple Power Analysis (SPA) and timing attacks. It evaluates exactly one unified addition and one dedicated doubling per scalar bit.
4. **Fast Cofactor Clearing (Receiver)**: Instead of multiplying by the massive full cofactor $h$, the receiver utilizes the drastically shorter scalar $\vert{}1-z\vert{}$ (approx. 80–90 bits).
5. **Hardware Parallelization Readiness**: The structural decoupling of the integer-based scalar recoding and field-based lookup tables exposes multiple vectors for Thread-Level and Instruction-Level Parallelism in compiled environments.

## Empirical Performance (Magma V2.29-6)

Benchmarking 100,000 deterministic iterations reveals that our strictly projective, inversion-free Jacobi pipeline completely absorbs the field-scaling penalty. Remarkably, our *interpreted* BLS12-479+ implementation in High-Throughput mode outperforms Magma's *natively compiled* baseline for BLS12-381 by over **20.2%**. The Montgomery ladder ensures constant-time execution while maintaining high efficiency:

| Curve Variant | Mode | Obfuscation (Cycles) | Deobfuscation (Cycles) |
| :--- | :--- | :--- | :--- |
| **BLS12-381** (Native) | Variable-Time | 26,432,936 | 2,666,708 |
| **BLS12-479+** (Interpreted) | **High-Throughput ($w\text{NAF}$)** | **21,088,209** | **2,632,178** |
| **BLS12-539+** (Interpreted) | **High-Throughput ($w\text{NAF}$)** | **25,758,478** | **3,060,789** |
| **BLS12-479+** (Interpreted) | **SPA-Resistant (Mont. Ladder)** | 34,516,524 | 5,782,380 |
| **BLS12-539+** (Interpreted) | **SPA-Resistant (Mont. Ladder)** | 40,939,151 | 6,084,008 |

*(Benchmarks executed on a 2.7 GHz Quad-Core Intel Core i7, 16GB RAM)*

## Repository Structure

* `bls12_479_g1_jacobi.m`: The full end-to-end obfuscation and deobfuscation pipeline test for the 479-bit target curve. Includes runtime verification of the logic and base clock-cycle benchmarking.
* `bls12_539_g1_jacobi.m`: The identical pipeline adapted for the ultra-high-security 539-bit target curve.
* `benchmarks++.m`: A standalone performance script evaluating 100,000 deterministic iterations of the High-Throughput $w\text{NAF}$ profile.
* `montgomery_ledder.m`: A standalone performance script evaluating 100,000 deterministic iterations of the SPA-Resistant Constant-Time Montgomery ladder profile.

## Usage

To run the scripts, you need the Magma Computer Algebra System installed. Run the scripts via the command line interface:

```bash
# Run the standalone benchmarks
magma benchmarks++.m
magma montgomery_ledder.m

# Run the detailed pipeline tests
magma bls12_479_g1_jacobi.m
magma bls12_539_g1_jacobi.m
```

### 📚 Academic Citation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23536594.svg)](https://doi.org/10.5281/zenodo.23536594)

If you utilize these scripts or the BLS12-479+ curve parameters in your research,
please cite the associated archive:

```bash
@misc{chmora2026bls12steganography,
  author       = {Andrey Chmora},
  title        = {Optimizing Elligator Squared for BLS Curves: A Jacobi Quartic Approach},
  month        = {Oct},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = 10.5281/zenodo.23536594
  url          = [https://zenodo.org/badge/DOI/10.5281/zenodo.23536594.svg](https://doi.org/10.5281/zenodo.23536594)
```
