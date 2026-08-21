#include "tise.hpp"

#define _USE_MATH_DEFINES

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <vector>
#include <map>
#include <regex>
#include "muParser.h"
#include <sstream>

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

std::vector<Real> buildEnergyGrid(Real E_threshold, Real E_max, int N_E)
{
    std::vector<Real> grid(N_E);
    const Real span = E_max - E_threshold;
    for (int i = 1; i <= N_E; ++i)
        grid[i - 1] = E_threshold + span * static_cast<Real>(i) / static_cast<Real>(N_E);
    return grid;
}

double radialPotential(double x, int L)
{
    const double l = static_cast<double>(L);
    return l * (l + 1.0) / (2.0 * x * x) - 1.0 / x;
}

// given a string, e.g. [0, 20), return true if input x falls in the bounds and false otherwise
// Bounds may be finite numbers or (+/-)inf/infinity (case-insensitive); brackets
// select inclusive ([, ]) vs. exclusive ((, )) endpoints.
bool inInterval(double x, const std::string& interval)
{
    static const std::regex re(
        R"(^\s*([\[\(])\s*(-?(?:\d+(?:\.\d+)?)|[+-]?(?:inf|infinity))\s*,\s*(-?(?:\d+(?:\.\d+)?)|[+-]?(?:inf|infinity))\s*([\]\)])\s*$)",
        std::regex_constants::icase
    );

    std::smatch m;
    if (!std::regex_match(interval, m, re)) {
        throw std::runtime_error("Invalid interval: " + interval);
    }

    auto parseBound = [](const std::string& s) {
        std::string t;
        t.reserve(s.size());
        for (char c : s)
            t += std::tolower(static_cast<unsigned char>(c));

        if (t == "inf" || t == "+inf" ||
            t == "infinity" || t == "+infinity")
            return std::numeric_limits<double>::infinity();

        if (t == "-inf" || t == "-infinity")
            return -std::numeric_limits<double>::infinity();

        return std::stod(s);
    };

    bool leftInclusive = m[1] == "[";
    bool rightInclusive = m[4] == "]";

    double left = parseBound(m[2]);
    double right = parseBound(m[3]);

    bool leftOK = leftInclusive ? (x >= left) : (x > left);
    bool rightOK = rightInclusive ? (x <= right) : (x < right);

    return leftOK && rightOK;
}

double evaluateFunction(std::map<std::string, std::string> function, double x)
{
    // given a map from domain to function, evaluate function for input x
    // for now, this naively assumes that the first match is the correct one
    // TODO: add error checking in case input x fits in multiple pieces, or function is not defined in terms of x
    for (const auto& [domain, fn] : function) {
        if (inInterval(x, domain))
        {
            // evaluate function
            // muparser expression using "x" as the sole variable, bound to a
            // fresh Parser each call since the expression string varies per piece
            mu::Parser p;
            p.DefineVar("x", &x);
            p.SetExpr(fn);

            // return result
            return p.Eval();
        }
    }
    std::ostringstream oss;
    oss << "Function domain does not cover x = "
        << x ;
    throw std::runtime_error(oss.str());
}

std::pair<std::vector<Real>, std::vector<Real>>
fillBandedMatrices(const bspline::BSpline &bs, int nEn, int order, int L, std::map<std::string, std::string> potential)
{
    std::vector<Real> Hmat(order * (nEn + 1), 0.0);
    std::vector<Real> Smat(order * (nEn + 1), 0.0);

    bspline::D2DFun fUni = [](double, const double *) { return 1.0; };
    // Piecewise potential supplied by the caller, evaluated per-x via muparser.
    // `L` is no longer used to select the potential here; it is retained for
    // eigenvalueError()'s comparison against the analytic hydrogen spectrum.
    bspline::D2DFun fPot = [potential](double x, const double *) {
        return evaluateFunction(potential, x);
    };
    // Previous hardcoded radial hydrogen-like potential, kept for reference:
    // bspline::D2DFun fPot = [L](double x, const double *) {
    //     return radialPotential(x, L);
    // };
    double parvec[1] = {0.0};

    auto bandIndex = [&](int row, int col) {
        return (row - 1) + (col - 1) * order;
    };

    for (int iBs2 = 2; iBs2 <= nEn + 2; ++iBs2)
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

std::pair<std::vector<Real>, std::vector<Real>> precomputeBoundaryCoupling(
    int order, int nEn, 
    std::vector<Real> Hmat, 
    std::vector<Real> Smat, 
    EigenResult eigen)
{
    // Calculate <phi_n | H | B_N> and <phi_n | B_N>. 
    // Note:
    // - H is in the basis made up by the B-Splines, hence H |B_N> is simply the N-th column of H.
    // - Hmat is not an N X N flattened matrix; it stores bands of width ``order`` from an ``nEn``x``nEn`` matrix, so there are ``order * nEn`` elements
    // - Hmat was made in column-major order (see fillBandedMatrices)
    // - eigen.vectors is a flattened set of all eigenvectors in "column-major" order
    // Thus we get H | B_N> by getting the last column of H.
    // Also note:
    // - Eigenvectors are NOT in the B-Spline basis, so we can't just compute <phi_i|B_N> by taking the last element of each |phi_n>.
    // instead we need to use SB_N, the last column of S, which accounts for overlap of all B-Splines with B_N.

    std::vector<Real> HB_N(nEn), SB_N(nEn);
    int iBs2 = nEn + 2;      // last B-spline index
    int col  = iBs2 - 1;     // = nEn + 1, matches bandIndex's 1-indexed col

    int iBs1Min = std::max(2, iBs2 - order + 1);
    for (int iBs1 = iBs1Min; iBs1 <= iBs2 - 1; ++iBs1) {
        int row = iBs1 + order - iBs2;            // band-local row, matches fill loop
        int idx = (row - 1) + (col - 1) * order;  // bandIndex(row, col)
        HB_N[iBs1 - 2] = Hmat[idx];               // iBs1=2 -> index 0
        SB_N[iBs1 - 2] = Smat[idx];               // same band position, from Smat
    }


    // <phi_n|H|B_N> / <phi_n|B_N> are the scalar product of phi_n with HB_N / SB_N.
    std::vector<Real> coeffs1(nEn), coeffs2(nEn);
    for (int i = 0; i < nEn; ++i)
    {
        Real hSum = 0.0, sSum = 0.0;
        for (int j = 0; j < nEn; ++j)
        {
            Real c = eigen.vectors[i * eigen.ldz + j];
            hSum += c * HB_N[j];
            sSum += c * SB_N[j];
        }
        coeffs1[i] = hSum;
        coeffs2[i] = sSum;
    }

    return {coeffs1, coeffs2};
}

std::vector<std::vector<Real>> buildContinuumState(
    int order, int nEn, 
    std::vector<Real> Hmat, 
    std::vector<Real> Smat, 
    EigenResult eigen,
    std::vector<Real> grid)
{
    auto [coeffs1, coeffs2] = precomputeBoundaryCoupling(order, nEn, Hmat, Smat, eigen);
    
    std::vector<std::vector<Real>> states(grid.size(), std::vector<Real>(eigen.values.size() + 1, 0.0));

    // Compute the coeffs for each point on the energy grid
    for(int E_idx = 0; E_idx < grid.size(); ++E_idx)
    {
        for(int i = 0; i < nEn; ++i)
        {
            states[E_idx][i] = (coeffs1[i] - (grid[E_idx] * coeffs2[i])) / (grid[E_idx] - eigen.values[i]);
        }
        states[E_idx][nEn] = 1.0;
    }

    return states;

}

AsymptoticResult matchAsymptotic(
    const bspline::BSpline &bs, 
    std::vector<std::vector<Real>> states, 
    const EigenResult &eigen, 
    std::vector<Real> grid, Real R
)
{
    AsymptoticResult result;
    result.A_E = std::vector<Real>(grid.size(), 0.0);
    result.delta = std::vector<Real>(grid.size(), 0.0);
    result.dDeltaDE = std::vector<Real>(grid.size(), 0.0);

    int nEn = eigen.dim;
    int nBSplines = bs.getNBSplines();

    // for each E, first find \bar psi_E(R) and \bar psi'_E(R), then calculate A_E, delta
    for (int E_idx = 0; E_idx < grid.size(); ++E_idx)
    {
        // states[E_idx][i] (i=0..nEn-1) give the coefficients for bound eigenstate phi_i which builds |\bar psi_E>, 
        // (and states[E_idx][nEn] is always 1.0, the coefficient for B_N)
        // To compute the scalar product needed to find \bar psi_E(R) and \bar psi'_E(R) with bs.eval, however, 
        // we need ALL coefficients to correspond to the B-Splines, not phi_n!!!
        std::vector<Real> fc(nBSplines, 0.0);
        for (int j = 0; j < nEn; ++j)
        {
            Real coeff = 0.0;
            for (int n = 0; n < nEn; ++n)
                coeff += states[E_idx][n] * eigen.vectors[n * eigen.ldz + j];
            fc[1 + j] = coeff;
        }
        fc[nEn + 1] = states[E_idx][nEn];

        // NOW we can use bs.eval
        Real psi_R = bs.eval(R, fc.data(), fc.size(), 0);
        Real psiPrime_R = bs.eval(R, fc.data(), fc.size(), 1);

        Real k = sqrt(2 * grid[E_idx]);

        result.A_E[E_idx] = sqrt(
            (2 / M_PI) / 
            (
                k * pow(psi_R, 2) +
                pow(psiPrime_R, 2) / k
            )
        );

        result.delta[E_idx] = std::atan(
            (k * psi_R) /
            (psiPrime_R)
        ) - (k * R);
    }

    // now that result.delta is filled, we can find result.dDeltaDE
    for (int E_idx = 0; E_idx < grid.size(); ++E_idx)
    {
        Real dE = grid[1] - grid[0];
        Real dSin2DeltaDE;
        if (E_idx == 0)
            dSin2DeltaDE = (std::sin(2*result.delta[E_idx+1]) - std::sin(2*result.delta[E_idx])) / dE;
        else if (E_idx == grid.size() - 1)
            dSin2DeltaDE = (std::sin(2*result.delta[E_idx]) - std::sin(2*result.delta[E_idx-1])) / dE; // TODO: will this sign be wrong
        else
            dSin2DeltaDE = (std::sin(2*result.delta[E_idx+1]) - std::sin(2*result.delta[E_idx-1])) / (2*dE);

        result.dDeltaDE[E_idx] = dSin2DeltaDE / (2.0 * std::cos(2.0 * result.delta[E_idx]));
    }

    return result;
}

void writeContinuumInfo(std::ostream &out,
                     const bspline::BSpline &bs,
                     const AsymptoticResult &result,
                     const std::vector<Real> &grid,
                     const std::vector<std::vector<Real>> &states,
                     std::vector<std::ostream *> stateOut,
                     int npts,
                     Real rMin,
                     Real rMax)
{
    // phase_shifts.dat: 3-col epsilon_i, delta(epsilon_i), dDelta/dE
    out << std::scientific << std::setprecision(16);
    for (std::size_t i = 0; i < grid.size(); ++i)
    {
        out << " " << std::setw(24) << grid[i]
            << " " << std::setw(24) << result.delta[i]
            << " " << std::setw(24) << result.dDeltaDE[i] << "\n";
    }

    // continuum_state_NNN.dat: 2-col x, psi_{epsilon_i}(x), one block per energy
    for (std::size_t i = 0; i < grid.size(); ++i)
    {
        std::ostream &stOut = *stateOut[i];
        stOut << std::scientific << std::setprecision(16);
        const int n = static_cast<int>(states[i].size());
        for (int ix = 1; ix <= npts; ++ix)
        {
            Real x = rMin + (rMax - rMin) *
                     static_cast<Real>(ix - 1) / static_cast<Real>(npts - 1);
            Real psi = result.A_E[i] * bs.eval(x, states[i].data(), n);
            stOut << " " << std::setw(24) << x
                  << " " << std::setw(24) << psi << "\n";
        }
    }
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

EigenResult solveTISE(int nNodes, int order, Real rMin, Real rMax, int L, std::map<std::string, std::string> potential,
                       Real E_threshold, Real E_max, int N_E)
{
    auto grid = buildUniformRadialGrid(nNodes, rMin, rMax);

    bspline::BSpline bs;
    int initInfo = bs.init(nNodes, order, grid);
    if (initInfo != 0)
        throw std::runtime_error("BSpline::init failed with code " +
                                 std::to_string(initInfo));

    int nBSplines = bs.getNBSplines();
    int nEn       = nBSplines - 2;

    auto [H, S] = fillBandedMatrices(bs, nEn, order, L, potential);
    EigenResult er = solveGeneralizedEigenproblem(H, S, nEn, order);

    auto energyGrid = buildEnergyGrid(E_threshold, E_max, N_E);
    std::vector<std::vector<tise::Real>> states = buildContinuumState(order, nEn, H, S, er, energyGrid);
    AsymptoticResult ar = matchAsymptotic(bs, states, er, energyGrid, 100);

    std::ofstream phaseShiftsOut("phase_shifts.dat");
    std::vector<std::ofstream> continuumStateFiles;
    continuumStateFiles.reserve(energyGrid.size());
    std::vector<std::ostream *> continuumStateOut;
    for (std::size_t i = 0; i < energyGrid.size(); ++i)
    {
        std::ostringstream oss;
        oss << "continuum_state_" << std::setw(3) << std::setfill('0') << (i + 1) << ".dat";
        continuumStateFiles.emplace_back(oss.str());
        continuumStateOut.push_back(&continuumStateFiles.back());
    }
    writeContinuumInfo(phaseShiftsOut, bs, ar, energyGrid, states, continuumStateOut, 301, rMin, rMax);

    return er;
}

} // namespace tise
