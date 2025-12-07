// BSpline C++ module
//
// This file is a C++ translation and adaptation of the original Fortran
// ModuleBSpline code by Luca Argenti.
//
// ---------------------------------------------------------------------
// Original copyright and license notice:
//
// ModuleBSplines, Copyright (C) 2020 Luca Argenti, PhD - Some Rights Reserved
// ModuleBSplines is licensed under a
// Creative Commons Attribution-ShareAlike 4.0 International License.
// A copy of the license is available at <http://creativecommons.org/licenses/by-nd/4.0/>.
//
// Luca Argenti is Associate Professor of Physics, Optics and Photonics
// Department of Physics and the College of Optics
// University of Central Florida
// 4111 Libra Drive, Orlando, Florida 32816, USA
// email: luca.argenti@ucf.edu
//
// ---------------------------------------------------------------------
//
// C++ translation and modifications:
//   Copyright (C) 2025  Martin de Salterain
//   Copyright (C) 2025  Hui Xian Grace Lim
//
// This C++ version is a derivative work based on ModuleBSplines and is
// distributed under the same Creative Commons license as the original.
// Please refer to the original license terms for details.

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include <iomanip> // for std::setw, std::setprecision, std::scientific

#include "BSpline.hpp"

// ------------------------------------------------------------------
// Select angular momentum quantum number l.
// The radial potential is taken as:
//
//   V_l(x) = l (l + 1) / (2 x^2) - 1 / x
//
// and the exact hydrogenic energies are
//
//   E_n = -1 / [ 2 (n_eff)^2 ],   where  n_eff = iEn + l.
//
// So:
//   l = 0:  V = -1/x,            E_n = -1 / [2 (iEn)^2]
//   l = 1:  V = 1/x^2 - 1/x,     E_n = -1 / [2 (iEn+1)^2]
// ------------------------------------------------------------------
constexpr int L = 1; // <-- change this integer to 0, 1, 2, ...

// Constants for maths...
const double hbar = 1.0; // in atomic units
const double pi = 3.14159265358979;
const double m_e = 3.75; // in atomic units
const double initial_position = 10.0; // initial position of the wavepacket

// Fortran LAPACK routine DSBGV
extern "C"
{
     void dsbgv_(char *jobz, char *uplo,
                 int *n, int *ka, int *kb,
                 double *ab, int *ldab,
                 double *bb, int *ldbb,
                 double *w,
                 double *z, int *ldz,
                 double *work,
                 int *info);
}

// Unity operator used in the problem (Template.f90)
double Unity(double /*x*/, const double * /*parvec*/)
{
     return 1.0;
}

// Radial potential for angular momentum l:
//
//   V_l(x) = l(l+1) / (2 x^2) - 1/x.
//
// For x > 0 (radial coordinate), and integer L >= 0.
double Potential(double x, const double * /*parvec*/)
{
     const double l = static_cast<double>(L);
     return l * (l + 1.0) / (2.0 * x * x) - 1.0 / x;
}

int main()
{
     // Bring a few types into local scope for convenience:
     //  - BSpline : the B-spline basis class
     //  - Real    : typedef for double
     //  - D2DFun  : std::function<double(double,const double*)> for local operators
     using bspline::BSpline;
     using bspline::D2DFun;
     using bspline::Real;

     // ------------------------------------------------------------
     // 1. Define B-spline grid and basis parameters
     // ------------------------------------------------------------
     // BS_NNODS : number of *distinct* radial grid points (nodes) in [BS_GRMIN, BS_GRMAX]
     // BS_ORDER : B-spline order k. Each B-spline is a piecewise polynomial of degree k-1
     //            that spans k intervals. Here k = 12 (order-12 splines).
     const int BS_NNODS = 51;
     const int BS_ORDER = 12;

     // Radial domain: [BS_GRMIN, BS_GRMAX]
     // This is the box in which we solve the radial Schrödinger equation.
     const Real BS_GRMIN = 0.0;
     const Real BS_GRMAX = 100.0;
     const Real BS_INTER = BS_GRMAX - BS_GRMIN;

     // Construct a *uniform* grid of BS_NNODS points between BS_GRMIN and BS_GRMAX.
     // These are the *distinct* nodes; the BSpline class will internally extend this
     // with repeated end knots to implement the usual B-spline boundary conditions.
     std::vector<Real> BS_GRID(BS_NNODS);
     for (int iNode = 1; iNode <= BS_NNODS; ++iNode)
     {
          BS_GRID[iNode - 1] =
              BS_GRMIN + BS_INTER * static_cast<Real>(iNode - 1) /
                             static_cast<Real>(BS_NNODS - 1);
     }

     // ------------------------------------------------------------
     // 2. Initialize the B-spline basis
     // ------------------------------------------------------------
     BSpline bspline;
     int info = bspline.init(BS_NNODS, BS_ORDER, BS_GRID);
     if (info != 0)
     {
          // If initialization fails, print an error and exit.
          // Possible causes: invalid node count, non-monotonic grid, etc.
          std::cerr << "Error in BSpline::init, info = " << info << "\n";
          return EXIT_FAILURE;
     }

     // Total number of B-splines is Nnodes + order - 2 (as in the Fortran code).
     const int nBSplines = bspline.getNBSplines(); // = BS_NNODS + BS_ORDER - 2

     // We impose that the radial wavefunction vanishes at the edges of the interval.
     // The first and last B-splines do not vanish at the boundaries, so we exclude them.
     // Therefore the variational subspace uses nEn = nBSplines - 2 B-splines.
     const int nEn = nBSplines - 2;                // exclude first and last B-spline

     // ------------------------------------------------------------
     // 3. Allocate banded Hamiltonian and overlap matrices
     // ------------------------------------------------------------
     // We store both H and S in *symmetric banded* storage, as required by LAPACK's DSBGV.
     //
     // Layout for a symmetric band matrix A of size n with upper bandwidth (ka):
     //   - We use a 2D array AB(ka+1, n) in column-major form.
     //   - The diagonal A(i,i) is stored at AB(ka+1, i).
     //   - The k-th superdiagonal A(i, i+k) is stored at AB(ka+1-k, i+k).
     //
     // Here BS_ORDER - 1 is the half-bandwidth ka = kb.
     //
     // We allocate Hmat and Smat as 1D vectors of length BS_ORDER * nEn,
     // which we will treat as (BS_ORDER x nEn) column-major arrays.
     std::vector<Real> Hmat(BS_ORDER * nEn, 0.0); // Hamiltonian
     std::vector<Real> Smat(BS_ORDER * nEn, 0.0); // Overlap (mass) matrix

     // ------------------------------------------------------------
     // 4. Fill H and S matrices
     // ------------------------------------------------------------
     {
          // Local operator functors:
          //  - fUni(x)  = 1            (used for overlap & kinetic energy kernel)
          //  - fPot(x)  = V_l(x)       (radial potential)
          D2DFun fUni = Unity;
          D2DFun fPot = Potential;

          // The Fortran interface for potential allowed a parvec array; we keep the
          // same interface for compatibility, even though we do not use it here.
          double parvec[1] = {0.0};

          // Small helper that maps (row, col) in 1-based Fortran-like band storage
          // into a 0-based index into our 1D vector.
          auto bandIndex = [](int row, int col, int ld)
          {
               // row : 1..ld   (ld = BS_ORDER)
               // col : 1..nEn  (matrix column index)
               // ld  : leading dimension (number of rows = BS_ORDER)
               return (row - 1) + (col - 1) * ld;
          };

          // Loop over column index iBs2 in the reduced basis (2 .. nEn+1).
          // In the full B-spline set, indices run from 1 .. nBSplines.
          // Our working basis uses indices 2 .. nBSplines-1, so the
          // generalized eigenproblem is of size nEn.
          for (int iBs2 = 2; iBs2 <= nEn + 1; ++iBs2)
          {
               // For a given column (B-spline index iBs2), only rows (iBs1)
               // that are within (BS_ORDER - 1) of iBs2 contribute to the band,
               // i.e., |iBs1 - iBs2| < BS_ORDER.
               int iBs1Min = std::max(2, iBs2 - BS_ORDER + 1);
               for (int iBs1 = iBs1Min; iBs1 <= iBs2; ++iBs1)
               {
                    // Overlap matrix element S_{iBs1,iBs2} = <B_iBs1 | B_iBs2>
                    Real overlap = bspline.integral(fUni, iBs1, iBs2);

                    // Kinetic energy: T = -1/2 d^2/dr^2
                    // The matrix element is ∫ B'_i Bs'_j dr / 2 (with appropriate radial measure).
                    // Here bspline.integral(fUni, iBs1, iBs2, 1, 1) computes
                    // ∫ B_i^(1)(r) B_j^(1)(r) dr (with the internal r^2 factor in the module);
                    // dividing by 2 gives the -1/2 factor in the Hamiltonian.
                    Real kinetic = bspline.integral(fUni, iBs1, iBs2,
                                                    /*n1=*/1, /*n2=*/1) /
                                   2.0;

                    // Potential energy: V_l(r) = l(l+1)/(2 r^2) - 1/r
                    // The matrix element is ∫ B_i(r) V_l(r) B_j(r) dr (again with r^2 from the module).
                    Real potential = bspline.integral(fPot, iBs1, iBs2,
                                                      /*n1=*/0, /*n2=*/0,
                                                      parvec);

                    // Convert the pair of B-spline indices (iBs1,iBs2) to
                    // the band-storage coordinates (row, col) of the reduced nEn x nEn matrix:
                    //   col = iBs2 - 1  -> column index in 1..nEn
                    //   row = iBs1 + BS_ORDER - iBs2  -> row index in 1..BS_ORDER
                    int col = iBs2 - 1;               // 1..nEn
                    int row = iBs1 + BS_ORDER - iBs2; // 1..BS_ORDER
                    int idx = bandIndex(row, col, BS_ORDER);

                    // Store overlap and Hamiltonian entries.
                    Smat[idx] = overlap;
                    Hmat[idx] = kinetic + potential;
               }
          }
     }

     // ------------------------------------------------------------
     // 5. Solve the generalized eigenproblem H c = E S c using DSBGV
     // ------------------------------------------------------------
     {
          // Threshold for deciding when an eigenvalue is "accurate"
          // relative to the analytical hydrogen-like formula.
          const Real ERROR_THRESHOLD = 1.0e-10;

          // Number of radial points at which to sample each eigenfunction
          // when writing it to a file.
          const int npts = 301;

          // LAPACK will overwrite A and B, so we work on copies of H and S.
          std::vector<Real> A = Hmat;
          std::vector<Real> B = Smat;

          // eval: eigenvalues E_i
          std::vector<Real> eval(nEn, 0.0);       // eigenvalues

          // evec: eigenvectors Z (dimension nEn x nEn, column-major).
          //       Column iEn corresponds to the iEn-th eigenvector.
          int ldz = nEn;                          // no extra padding in C++
          std::vector<Real> evec(ldz * nEn, 0.0); // eigenvectors

          // Workspace array of length 3*nEn as required by DSBGV.
          std::vector<Real> work(3 * nEn, 0.0);

          // LAPACK DSBGV arguments:
          //   jobz = 'V' -> compute eigenvalues and eigenvectors
          //   uplo = 'U' -> upper triangle of band matrices is stored
          char jobz = 'V';
          char uplo = 'U';
          int n = nEn;           // dimension of the generalized eigenproblem
          int ka = BS_ORDER - 1; // number of superdiagonals in H
          int kb = BS_ORDER - 1; // number of superdiagonals in S
          int ldab = BS_ORDER;   // leading dimension of A (rows in band storage)
          int ldbb = BS_ORDER;   // leading dimension of B
          int infoLapack = 0;

          // Call LAPACK DSBGV:
          //
          //   A and B are symmetric band matrices (upper storage).
          //   On exit, eval contains eigenvalues in ascending order,
          //   and evec contains the corresponding eigenvectors (columns of Z).
          dsbgv_(&jobz, &uplo,
                 &n, &ka, &kb,
                 A.data(), &ldab,
                 B.data(), &ldbb,
                 eval.data(),
                 evec.data(), &ldz,
                 work.data(),
                 &infoLapack);

               int kd = BS_ORDER - 1;
               int ld = BS_ORDER;

               if (infoLapack != 0)
               {
                    std::cerr << "DSBGV info: " << infoLapack << "\n";
                    return EXIT_FAILURE;
               }

          // Print numeric values in scientific notation with high precision.
          std::cout << std::scientific << std::setprecision(16);

          int iEn;
          for (iEn = 1; iEn <= nEn; ++iEn)
          {
               // Eigenvalue for state iEn (LAPACK gives them in increasing order).
               Real eig = eval[iEn - 1];

               // ------------------------------------------------------------
               // 5a. Compute analytic hydrogenic reference energy and error
               // ------------------------------------------------------------
               // Effective principal quantum number:
               //   n_eff = iEn + L
               // This matches the logic:
               //   L = 0: n_eff = iEn
               //   L = 1: n_eff = iEn + 1
               Real n_eff = static_cast<Real>(iEn + L);

               // Analytic hydrogenic energy: E_exact = -1 / [2 n_eff^2].
               // Here we store the *difference* eig - E_exact as:
               //   eigenvalueError = eig + 1/(2 n_eff^2)
               // (so small |eigenvalueError| means good agreement).
               Real eigenvalueError = eig + 1.0 / (2.0 * n_eff * n_eff);

               // If the error is larger than the threshold, we stop counting
               // further states as "accurate" (as the Fortran code does).
               if (eigenvalueError > ERROR_THRESHOLD)
               {
                    break;
               }

               // Print the eigenvalue and its deviation from the analytic value.
               std::cout << std::setw(4) << iEn << "  "
                         << std::setw(24) << eig << "  "
                         << std::setw(24) << eigenvalueError << "\n";

               // ------------------------------------------------------------
               // 5b. Write the eigenfunction ψ_iEn(r) to a file
               // ------------------------------------------------------------
               // File name: "EigenState_XXX" (e.g., EigenState_001, EigenState_002, ...)
               std::ostringstream oss;
               oss << "EigenState_" << std::setw(3) << std::setfill('0') << iEn;
               std::string filename = oss.str();

               std::ofstream out(filename);
               if (!out)
               {
                    std::cerr << "Cannot open file " << filename << " for writing\n";
                    return EXIT_FAILURE;
               }
               out << std::scientific << std::setprecision(16);

               // Each column of evec is a generalized eigenvector in the reduced B-spline basis
               // (indices 2 .. nBSplines-1). We extract column iEn as zCol.
               const double *zCol = &evec[(iEn - 1) * ldz];

               // Construct a full coefficient vector of length nBSplines for the BSpline::eval call:
               //   coeffs[0]          = 0       -> coefficient of B-spline #1
               //   coeffs[1..nEn]     = zCol[]  -> coefficients of B-splines #2 .. #nBSplines-1
               //   coeffs[nEn+1]      = 0       -> coefficient of B-spline #nBSplines
               //
               // This matches the Fortran code, where the first and last B-splines are not used
               // in the generalized eigenproblem but still exist in the BSpline set.
               std::vector<Real> coeffs(nBSplines, 0.0);
               for (int k = 0; k < nEn; ++k)
               {
                    coeffs[1 + k] = zCol[k];
               }

               // Sample the eigenfunction on npts uniformly spaced r-values in [BS_GRMIN, BS_GRMAX].
               for (int ix = 1; ix <= npts; ++ix)
               {
                    Real x = BS_GRMIN +
                             (BS_GRMAX - BS_GRMIN) *
                                 static_cast<Real>(ix - 1) /
                                 static_cast<Real>(npts - 1);

                    // Evaluate ψ(r) = sum_Bs coeffs[Bs-1] * B_spline_Bs(r)
                    Real fun = bspline.eval(x, coeffs.data(), nBSplines);

                    out << " " << std::setw(24) << x
                        << " " << std::setw(24) << fun << "\n";
               }
               out.close();
          }

          // iEn is incremented one past the last accepted index, so iEn-1
          // is the number of eigenvalues that passed the accuracy threshold.
          std::cout << "Number of Accurate Eigenvalues : " << (iEn - 1) << "\n";
     }

     return 0;
}
