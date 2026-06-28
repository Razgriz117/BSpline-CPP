#pragma once

#include <iosfwd>
#include <vector>
#include "BSpline.hpp"

namespace tise
{

using Real = bspline::Real;

// Holds the output of solveGeneralizedEigenproblem.
struct EigenResult
{
    std::vector<Real> values;  // eigenvalues, ascending, size dim
    std::vector<Real> vectors; // eigenvectors, column-major, size dim*dim
    int dim;                   // matrix dimension (nEn)
    int ldz;                   // leading dimension (== dim)
};

// Build a uniform radial grid of nNodes points on [rMin, rMax].
std::vector<Real> buildUniformRadialGrid(int nNodes, Real rMin, Real rMax);

// Radial hydrogen-like potential V_L(x) = L(L+1)/(2x^2) - 1/x.
double radialPotential(double x, int L);

// Fill symmetric banded Hamiltonian H and overlap S matrices (LAPACK 'U' storage).
// Returns {Hmat, Smat}, each of size order * nEn.
std::pair<std::vector<Real>, std::vector<Real>>
fillBandedMatrices(const bspline::BSpline &bs, int nEn, int order, int L);

// Solve H c = E S c via LAPACK DSBGV.
// H and S are consumed (overwritten); pass by value intentionally.
EigenResult solveGeneralizedEigenproblem(std::vector<Real> H,
                                         std::vector<Real> S,
                                         int nEn,
                                         int order);

// Analytic hydrogenic energy: E = -1 / (2 * (n + L)^2).
Real analyticHydrogenEnergy(int n, int L);

// Difference between computed eigenvalue and analytic energy for state n.
Real eigenvalueError(Real computed, int n, int L);

// Extract eigenvector column iEn (1-based) from the column-major evec array
// and embed it into a zero-padded vector of length nBSplines.
// coeffs[0] = 0, coeffs[1..nEn] = evec column, coeffs[nEn+1] = 0.
std::vector<Real> eigenstateCoefficients(const std::vector<Real> &evec,
                                          int iEn,
                                          int nEn,
                                          int nBSplines);

// Write a single eigenstate to `out`: npts lines of "x  psi(x)".
void writeEigenstate(std::ostream &out,
                     const bspline::BSpline &bs,
                     const std::vector<Real> &coeffs,
                     int npts,
                     Real rMin,
                     Real rMax);

// Top-level TISE solver: build grid, fill matrices, diagonalise.
EigenResult solveTISE(int nNodes, int order, Real rMin, Real rMax, int L);

} // namespace tise
