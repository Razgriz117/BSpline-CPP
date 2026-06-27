#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "utils.hpp"
#include "BSpline.hpp"

// Print a symmetric band matrix (LAPACK 'U' storage) as a dense n x n matrix.
//
// band  : band-storage array (column-major), length >= ld * n
// n     : matrix size (n x n)
// kd    : number of superdiagonals (half-bandwidth)
// ld    : leading dimension of the band matrix (number of rows = kd + 1, here BS_ORDER)
// name  : label to print before the matrix
// width : field width for each entry
// precision : number of digits after decimal (for scientific format)
void printBandMatrixAsDense(const std::vector<bspline::Real> &band,
                            int n,
                            int kd,
                            int ld,
                            const std::string &name,
                            int width,
                            int precision)
{
     using bspline::Real;
     using std::cout;
     using std::endl;

     cout << "\nMatrix " << name << " (dense " << n << "x" << n << " view):\n\n";

     // Configure numeric format: scientific with fixed width per entry.
     cout << std::scientific << std::setprecision(precision);

     // Helper lambda to get A(i,j) from band storage (1-based i,j).
     auto getA = [&](int i, int j) -> Real
     {
          // Exploit symmetry: always fetch from upper triangle (i <= j).
          if (i > j)
               std::swap(i, j);

          // If outside the band, the entry is zero.
          if (j - i > kd)
               return Real(0.0);

          // LAPACK formula for UPLO='U':
          //   AB(kd+1 + i - j, j) = A(i, j)
          int row = kd + 1 + i - j; // 1..kd+1
          int col = j;              // 1..n

          int idx = (row - 1) + (col - 1) * ld; // column-major indexing
          return band[idx];
     };

     // Loop over rows and columns of the dense matrix
     for (int i = 1; i <= n; ++i)
     {
          for (int j = 1; j <= n; ++j)
          {
               Real val = getA(i, j);
               cout << std::setw(width) << val;
          }
          cout << '\n';
     }

     cout << endl;
}

// Build C = Z^T from the eigenvector array `evec`.
//
// - evec: column-major matrix Z of size ldz x n (here ldz = n = nEn).
// - n   : number of eigenvectors / problem size.
// - ldz : leading dimension of Z (ldz >= n).
//
// Returns C as a column-major n x n matrix: C(i,j) = Z(j,i).
std::vector<bspline::Real> buildAdjointMatrix(const std::vector<bspline::Real> &evec,
                   int n,
                   int ldz)
{
     using bspline::Real;

     std::vector<Real> C(n * n, Real(0));

     // Z(i,j) = evec[i + j*ldz]  (0-based i,j)
     // C(i,j) = Z(j,i)           (0-based i,j)
     //      => C[i + j*n] = evec[j + i*ldz]
     for (int j = 0; j < n; ++j)
     { // column index of C
          for (int i = 0; i < n; ++i)
          { // row index of C
               C[i + j * n] = evec[j + i * ldz];
          }
     }

     return C;
}