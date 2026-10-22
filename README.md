# Optimizing Elligator Squared for BLS Curves via Jacobi Quartics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23537349.svg)](https://doi.org/10.5281/zenodo.23537349)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the Magma Computer Algebra System proof-of-concept implementation for the steganographic obfuscation of BLS12-family public keys. 

Specifically, this repository demonstrates a highly optimized Elligator Squared pipeline for the **BLS12-479+** and **BLS12-539+** curves. By leveraging their algebraically engineered 2-torsion subgroups, the pipeline constructs a mathematically exact birational equivalence to the **Extended Jacobi Quartic** space, effectively shattering the traditional computational penalty of public key obfuscation.

## Mathematical & Architectural Highlights

1. **Inversion-Free Arithmetic**: By carefully fusing the Billet-Joye birational mappings with Hişil's homogeneous projective coordinates $(X:Y:Z)$, we completely eliminate intermediate field inversions. The pipeline strictly defers the solitary field inversion to the final affine projection at the very end of the deobfuscation cycle.
2. **Jacobi Quartic Sandbox**: Heavy cryptographic operations are executed entirely within the accelerated Jacobi space, utilizing an optimized dedicated doubling formula ($3\mathbf{M}+4\mathbf{S}+1\mathbf{C}+3\mathbf{add}$).
3. **Asymmetric $w\text{NAF}$ with Batch-Normalized Table Generation**: To strictly optimize latency across scalars of vastly different bit lengths, the pipeline employs a wide window ($w=5$) for heavy full-scalar blinding and a minimal window ($w=2$) to avoid table overhead during cofactor clearing. Precomputation tables are normalized ($Z=1$) using Montgomery Batch Inversion, enabling the exclusive use of dedicated addition formulas ($9\mathbf{M}+5\mathbf{S}+1\mathbf{C}+6\mathbf{add}$) during scalar multiplication loops.
4. **Fast Cofactor Clearing (Receiver)**: Instead of multiplying by the massive full cofactor $h$, the receiver utilizes the drastically shorter scalar $\vert{}1-z\vert{}$ (approx. 80–90 bits). The $w\text{NAF}$ multiplication is performed fully inside the Jacobi space, keeping the receiver's verification overhead to an absolute minimum.
5. **Hardware Parallelization Readiness**: The structural decoupling of the integer-based scalar recoding and field-based lookup tables, alongside independent coordinate transformations during explicit 2-isogeny evaluations, exposes multiple vectors for Thread-Level and Instruction-Level Parallelism in compiled environments.

## Empirical Performance (Magma V2.29-6)

Benchmarking 100,000 deterministic iterations reveals that our strictly projective, inversion-free Jacobi pipeline completely absorbs the field-scaling penalty. Remarkably, our *interpreted* BLS12-479+ implementation outperforms Magma's *natively compiled* baseline for BLS12-381 by over **20.2%**, and even the 539-bit variant outperforms the 381-bit native baseline by **2.5%**:

| Curve Variant | Coordinate Space | Obfuscation (Cycles) | Deobfuscation (Cycles) |
| :--- | :--- | :--- | :--- |
| **BLS12-381** (Native) | Weierstrass (Mixed) | 26,432,936 | 2,666,708 |
| **BLS12-479+** (Interpreted) | **Jacobi Quartic ($\mathcal{P}$)** | **21,088,209** | **2,632,178** |
| **BLS12-539+** (Interpreted) | **Jacobi Quartic ($\mathcal{P}$)** | **25,758,478** | **3,060,789** |

*(Benchmarks executed on a 2.7 GHz Quad-Core Intel Core i7, 16GB RAM)*

## Repository Structure

* `bls12_479_g1_jacobi.m`: The full end-to-end obfuscation and deobfuscation pipeline test for the 479-bit target curve. Includes runtime verification of the logic and base clock-cycle benchmarking.
* `bls12_539_g1_jacobi.m`: The identical pipeline adapted for the ultra-high-security 539-bit target curve.
* `benchmarks.m`: A standalone performance script designed to run headless. It evaluates 100,000 deterministic iterations for both curves to output clean CPU clock-cycle metrics.

## Usage

To run the scripts, you need the Magma Computer Algebra System installed. Run the scripts via the command line interface:

```bash
# Run the standalone benchmarks for both curves
magma benchmarks.m

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
