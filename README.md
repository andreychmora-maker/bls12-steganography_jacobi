# Optimizing Elligator Squared for BLS Curves via Jacobi Quartics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23537349.svg)](https://doi.org/10.5281/zenodo.23537349)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the Magma Computer Algebra System proof-of-concept implementation for the steganographic obfuscation of BLS12-family public keys. 

Specifically, this repository demonstrates a highly optimized Elligator Squared pipeline for the **BLS12-479+** and **BLS12-539+** curves. By leveraging their algebraically engineered 2-torsion subgroups, the pipeline constructs a mathematically exact birational equivalence to the **Extended Jacobi Quartic** space, effectively shattering the traditional computational penalty of public key obfuscation.

## Mathematical & Architectural Highlights

1. **Inversion-Free Arithmetic**: By carefully fusing the Billet-Joye birational mappings with Hişil's homogeneous projective coordinates $(X:Y:Z)$, we completely eliminate intermediate field inversions. The pipeline strictly defers the solitary field inversion to the final affine projection at the very end of the deobfuscation cycle.
2. **Jacobi Quartic Sandbox**: Heavy cryptographic operations are executed entirely within the accelerated Jacobi space, utilizing an optimized dedicated doubling formula ($3\mathbf{M}+4\mathbf{S}+1\mathbf{C}+3\mathbf{add}$).
3. **Fast Cofactor Blinding (Sender)**: Full-scalar blinding is executed inside the Jacobi space using $w\text{NAF}$ (Non-Adjacent Form) recoding.
4. **Fast Cofactor Clearing (Receiver)**: Instead of multiplying by the massive full cofactor $h$, the receiver utilizes the drastically shorter scalar $\vert{}1-z\vert{}$ (approx. 80-90 bits). The $w\text{NAF}$ multiplication is performed fully inside the Jacobi space, keeping the receiver's verification overhead to an absolute minimum.

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
  doi          = 10.5281/zenodo.23505994
  url          = [https://zenodo.org/badge/DOI/10.5281/zenodo.23505994.svg](https://doi.org/10.5281/zenodo.23505994)
```
