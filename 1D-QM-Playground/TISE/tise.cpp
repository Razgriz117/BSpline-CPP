#include "tise.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <vector>

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

namespace tise
{

std::vector<Real> buildUniformRadialGrid(int nNodes, Real rMin, Real rMax)
{
    std::vector<Real> grid(nNodes);
    const Real span = rMax - rMin;
    for (int i = 0; i < nNodes; ++i)
        grid[i] = rMin + span * static_cast<Real>(i) / static_cast<Real>(nNodes - 1);
    return grid;
}

double radialPotential(double x, int L)
{
    const double l = static_cast<double>(L);
    return l * (l + 1.0) / (2.0 * x * x) - 1.0 / x;
}

std::pair<std::vector<Real>, std::vector<Real>>
fillBandedMatrices(const bspline::BSpline &bs, int nEn, int order, int L)
{
    std::vector<Real> Hmat(order * nEn, 0.0);
    std::vector<Real> Smat(order * nEn, 0.0);

    bspline::D2DFun fUni = [](double, const double *) { return 1.0; };
    bspline::D2DFun fPot = [L](double x, const double *) {
        return radialPotential(x, L);
    };
    double parvec[1] = {0.0};

    auto bandIndex = [&](int row, int col) {
        return (row - 1) + (col - 1) * order;
    };

    for (int iBs2 = 2; iBs2 <= nEn + 1; ++iBs2)
    {
        int iBs1Min = std::max(2, iBs2 - order + 1);
        for (int iBs1 = iBs1Min; iBs1 <= iBs2; ++iBs1)
        {
            Real overlap  = bs.integral(fUni, iBs1, iBs2);
            Real kinetic  = bs.integral(fUni, iBs1, iBs2, 1, 1) / 2.0;
            Real potential = bs.integral(fPot, iBs1, iBs2, 0, 0, parvec);

            int col = iBs2 - 1;
            int row = iBs1 + order - iBs2;
            int idx = bandIndex(row, col);

            Smat[idx] = overlap;
            Hmat[idx] = kinetic + potential;
        }
    }

    return {Hmat, Smat};
}

EigenResult solveGeneralizedEigenproblem(std::vector<Real> H,
                                          std::vector<Real> S,
                                          int nEn,
                                          int order)
{
    EigenResult result;
    result.dim = nEn;
    result.ldz = nEn;
    result.values.assign(nEn, 0.0);
    result.vectors.assign(nEn * nEn, 0.0);

    std::vector<Real> work(3 * nEn, 0.0);

    char jobz = 'V';
    char uplo = 'U';
    int n    = nEn;
    int ka   = order - 1;
    int kb   = order - 1;
    int ldab = order;
    int ldbb = order;
    int ldz  = nEn;
    int info = 0;

    dsbgv_(&jobz, &uplo,
            &n, &ka, &kb,
            H.data(), &ldab,
            S.data(), &ldbb,
            result.values.data(),
            result.vectors.data(), &ldz,
            work.data(),
            &info);

    if (info != 0)
        throw std::runtime_error("DSBGV failed with info=" + std::to_string(info));

    return result;
}

Real analyticHydrogenEnergy(int n, int L)
{
    const double n_eff = static_cast<double>(n + L);
    return -1.0 / (2.0 * n_eff * n_eff);
}

Real eigenvalueError(Real computed, int n, int L)
{
    return computed - analyticHydrogenEnergy(n, L);
}

std::vector<Real> eigenstateCoefficients(const std::vector<Real> &evec,
                                          int iEn,
                                          int nEn,
                                          int nBSplines)
{
    std::vector<Real> coeffs(nBSplines, 0.0);
    const Real *col = &evec[(iEn - 1) * nEn];
    for (int k = 0; k < nEn; ++k)
        coeffs[1 + k] = col[k];
    return coeffs;
}

void writeEigenstate(std::ostream &out,
                     const bspline::BSpline &bs,
                     const std::vector<Real> &coeffs,
                     int npts,
                     Real rMin,
                     Real rMax)
{
    out << std::scientific << std::setprecision(16);
    const int n = static_cast<int>(coeffs.size());
    for (int ix = 1; ix <= npts; ++ix)
    {
        Real x = rMin + (rMax - rMin) *
                 static_cast<Real>(ix - 1) / static_cast<Real>(npts - 1);
        Real psi = bs.eval(x, coeffs.data(), n);
        out << " " << std::setw(24) << x
            << " " << std::setw(24) << psi << "\n";
    }
}

EigenResult solveTISE(int nNodes, int order, Real rMin, Real rMax, int L)
{
    auto grid = buildUniformRadialGrid(nNodes, rMin, rMax);

    bspline::BSpline bs;
    int initInfo = bs.init(nNodes, order, grid);
    if (initInfo != 0)
        throw std::runtime_error("BSpline::init failed with code " +
                                 std::to_string(initInfo));

    int nBSplines = bs.getNBSplines();
    int nEn       = nBSplines - 2;

    auto [H, S] = fillBandedMatrices(bs, nEn, order, L);
    return solveGeneralizedEigenproblem(std::move(H), std::move(S), nEn, order);
}

} // namespace tise
