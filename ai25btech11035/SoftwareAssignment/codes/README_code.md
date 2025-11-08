# Matrix Theory (AI1000) – Course Project  
**Author:** AI25BTECH11035 — SUJAL RAJANI  
**Title:** Power Block Orthogonal Iteration (PBOI) for Image Compression


##  Overview

This project implements **Truncated Singular Value Decomposition (SVD)** using the **Power Block Orthogonal Iteration (PBOI)** method to compress grayscale images efficiently.

Unlike full SVD or Randomized SVD, this approach uses repeated multiplication of \(A^T A\) with a block of vectors, followed by orthonormalization, to extract the **dominant singular subspace** of the image matrix.

The goal is to obtain a low-rank approximation:

\[
A_k = U_k \Sigma_k V_k^T
\]

where \(k\) is small compared to the image dimension, thus significantly reducing storage while retaining image quality.


##  Algorithm

**Power Block Orthogonal Iteration (PBOI) for SVD**

**Input:**  
Matrix \(A \in \mathbb{R}^{m 	imes n}\), rank \(k\), iterations \(q\)  

**Output:**  
\(U_k, \Sigma_k, V_k\) — truncated SVD for low-rank approximation


### Pseudocode

Initialize Q ← random orthonormal matrix of size (n × k)

for i = 1 to q:
    Z ← A^T * (A * Q)             // Apply (A^T A) on Q
    Q ← QR_orthonormalize(Z)      // Modified Gram-Schmidt (MGS)

B ← A * Q                          // Projected matrix (m × k)

[U_tilde, Σ, V_tilde] ← SVD(B)    // SVD on smaller matrix

U ← U_tilde
V ← Q * V_tilde                   // Recover right singular vectors

Return U_k, Σ_k, V_k
##  Intuition

- \(A^T A\) is a symmetric positive semidefinite matrix.
- Its top eigenvectors correspond to the **right singular vectors** \(V_k\) of \(A\).
- By repeatedly multiplying \(A^T A\) with a random block \(Q\), the iteration converges to the span of the dominant eigenvectors.
- Each iteration refines \(Q\) and maintains orthonormality using QR decomposition.

Thus, the PBOI algorithm efficiently finds the same subspace as SVD — but **without** computing a full decomposition.


##  Implementation Steps

1. **Read grayscale image** (PGM format) → form matrix \(A\)
2. **Initialize random orthonormal block** \(Q\) of size (n × k)
3. **Repeat power block iteration (q times):**
   - Compute \(Z = A^T (A Q)\)
   - Orthonormalize columns of \(Z\) (Modified Gram-Schmidt)
   - Update \(Q ← Z\)
4. **Form projected matrix:** \(B = A Q\)
5. **Perform SVD of B:**  
   \(B = \widetilde{U} \Sigma \widetilde{V}^T\)
6. **Recover full singular vectors:**
   \[
   U = \widetilde{U}, \quad V = Q \widetilde{V}
   \]
7. **Reconstruct image:**
   \[
   A_k = U_k \Sigma_k V_k^T
   \]
8. **Write compressed image** (PGM format)



##  Compilation

gcc -O3 -o compress codes/c_main/main.c -lm




##  Running the Program

./compress

When prompted, provide:

- Input PGM filename  
- Output PGM filename  
- Rank \(k\)


## Example

| Parameter | Meaning | Typical Value |
|------------|----------|----------------|
| `K` | Target rank | 5, 10, 20, 50, 100 |
| `ITER` | Power iterations | 25 |
| `A` | Input image matrix | m × n |
| `Q` | Orthonormal block | n × K |


##  Output

- `output.pgm` — compressed image using top-\(k\) singular values  
- (Optional) `singular_values.txt` — list of top singular values


## Advantages of Power Block Method

| Feature | Description |
|----------|--------------|
| Iterative & memory efficient | No need for full SVD |
| Block orthogonalization | Stable for multiple vectors |
| Handles large matrices | Only needs matrix-vector products |
| Captures dominant features | Accurate low-rank reconstruction |


##  Summary

The **Power Block Orthogonal Iteration** algorithm provides an efficient way to approximate the dominant singular subspace of a large matrix by using iterative multiplication and orthonormalization.  

It achieves compression by retaining only the top-\(k\) singular components, reducing the image size while maintaining visual fidelity — an elegant combination of **linear algebra and numerical stability**.
