# SaberX4: High-throughput Software Implementationof Saber Key Encapsulation Mechanism
Paper available at https://eprint.iacr.org/2019/1309.pdf

Sujoy Sinha Roy

Abstract: Saber is a module lattice-based CCA-secure key encapsulation mechanism (KEM) which has been shortlisted for the second round of NIST's Post Quantum Cryptography Standardization project. To attain simplicity and efficiency on constrained devices, the Saber algorithm is serial by construction. However, on high-end platforms, such as modern Intel processors with AVX2 instructions, Saber achieves limited speedup using vector processing instructions due to its serial nature.

In this paper we overcome the above-mentioned algorithmic bottleneck and propose a high-throughput software implementation of Saber, which we call `SaberX4', targeting modern Intel processors with AVX2 vector processing support. We apply the batching technique at the highest level of the implementation hierarchy and process four Saber KEM operations simultaneously in parallel using the AVX2 vector processing instructions. Our proof-of-concept software implementation of SaberX4 achieves nearly 1.5 times higher throughput at the cost of latency degradation within acceptable margins, compared to the AVX2-optimized non-batched implementation of Saber by its authors.

We anticipate that both latency and throughput of SaberX4 will improve in the future with improved computer architectures and more optimization efforts.

