# Image Compression Using Truncated SVD

This repository contains the full C implementation of **image compression using Truncated Singular Value Decomposition (SVD)** as part of the *Software Assignment* for **AI25BTECH11035 - Sujal Rajani**.



## Author
**Sujal Rajani**  
Course: Software Assignment (AI1000)  
Submission: November 2025  



##  Project Structure

SoftwareAssignment/
├── codes/                        # Contains all C source files
│   ├── c_libs/                   # Header and helper C libraries
│   ├── c_main/                   # Main C program files
│   ├── python_driver/            # Python visualization and I/O scripts
│   └── hybrid_c_python/          # For mixed C + Python implementations
│       ├── c_backend/
│       └── python_frontend/
│
├── figs/                         # Original and reconstructed images
│
├── tables/                       # LaTeX tables (Frobenius norms, errors)
│   └── table.tex
│
├── report.pdf                    # Final IEEE-formatted project report
├── README.md                     # Project overview (this file)
└── README_code.md                # Detailed explanation of implementation


##  Project Summary

This project demonstrates how **Truncated SVD** can be used to compress grayscale images by approximating the original image matrix using only the top-\(k\) singular values and vectors.  
The implementation is written entirely in **C**, with optional **Python support** for visualization.

### The workflow includes:
- Reading and writing binary PGM images  
- Computing top-\(k\) SVD using **Power Iteration** and **QR decomposition**  
- Reconstructing compressed images  
- Calculating **Frobenius norm** and **relative reconstruction error**  




##  Key Concepts
- **Truncated SVD**: \( A \approx U_k \Sigma_k V_k^T \)
- **Power Iteration**: Efficiently approximates top eigenvectors of \( A^T A \)
- **Modified Gram-Schmidt QR**: Ensures orthonormality and numerical stability
- **Frobenius Norm**: Measures compression quality  
  \( \|A - A_k\|_F / \|A\|_F \)


##  Deliverables
- Full C source code (`codes/`)  
- Reconstructed images (`figs/`)  
- Error analysis and results (`tables/table.tex`)  
- Final IEEE-style report (`report.pdf`)  
- This README documentation  


##  Learning Outcome
Through this project, the following concepts were implemented and analyzed:
- Matrix factorization using iterative methods  
- Numerical stability using QR decomposition  
- Performance of compression vs. quality using Frobenius error metrics  
- Understanding trade-offs between rank \(k\) and reconstruction accuracy
