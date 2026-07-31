#include "tise.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
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

// === A1: boundary-condition asymptote classifier (REQ-F-030) ===
// See docs/planning/engineer-a-plan-A1.md for the derivation of the
// thresholds/window sizes below.
ConvergenceFit classifySequenceConvergence(const std::vector<Real> &V, Real ratio)
{
    const int N = static_cast<int>(V.size());

    std::vector<Real> dV(N - 1);
    Real maxAbsV = std::abs(V[N - 1]);
    Real maxAbsDV = 0.0;
    for (int k = 0; k < N - 1; ++k)
    {
        dV[k] = V[k + 1] - V[k];
        maxAbsV = std::max(maxAbsV, std::abs(V[k]));
        maxAbsDV = std::max(maxAbsDV, std::abs(dV[k]));
    }

    ConvergenceFit fit{};
    const Real nan = std::numeric_limits<Real>::quiet_NaN();

    // Flatness pre-check: successive differences are all numerically zero.
    if (maxAbsDV <= 1e-10 + 1e-9 * maxAbsV)
    {
        fit.isDivergent = false;
        fit.isFlat = true;
        fit.fittedLimit = V[N - 1];
        fit.powerLawExponent = nan;
        return fit;
    }

    // Divergence check: compare difference magnitude in the back half of the
    // window vs. the middle. If it hasn't shrunk substantially, the sequence
    // isn't converging to a finite limit.
    const int mid = (N - 2) / 2;
    const Real shrinkFactor = std::abs(dV[N - 2]) / std::max(std::abs(dV[mid]), 1e-300);
    if (shrinkFactor >= 0.5)
    {
        fit.isDivergent = true;
        fit.isFlat = false;
        fit.fittedLimit = nan;
        fit.powerLawExponent = nan;
        return fit;
    }

    // Power-law fit via successive-ratio over a tail window, avoiding both
    // near-field transients and far-field floating-point noise.
    const int windowStart = std::max(0, N - 7);
    const int windowEnd = N - 3; // inclusive; dV[windowEnd + 1] must stay in range
    std::vector<Real> pEstimates;
    for (int k = windowStart; k <= windowEnd && k + 1 < N - 1; ++k)
        pEstimates.push_back(std::log(std::abs(dV[k]) / std::abs(dV[k + 1])) / std::log(ratio));

    std::sort(pEstimates.begin(), pEstimates.end());
    const Real pFit = pEstimates[pEstimates.size() / 2]; // median

    const Real ratioP = std::pow(ratio, -pFit);
    const Real fittedLimit = (V[N - 1] - ratioP * V[N - 2]) / (1.0 - ratioP);

    fit.isDivergent = false;
    fit.isFlat = false;
    fit.fittedLimit = fittedLimit;
    fit.powerLawExponent = pFit;
    return fit;
}

namespace
{
constexpr Real kPi = 3.14159265358979323846;
}

Real case3WindowFunction(Real x, Real R, Real delta, DomainSide side)
{
    // Signed distance beyond the boundary: d > 0 means x is outside the box
    // (beyond the wall), d < 0 means x is inside the trusted region.
    const Real d = (side == DomainSide::Right) ? (x - R) : (R - x);

    if (d <= -delta)
        return 1.0;
    if (d >= 0.0)
        return 0.0;

    const Real theta = (kPi / 2.0) * (d / delta);
    const Real s = std::sin(theta);
    return s * s;
}

Real evaluateWindowedPotential(const std::map<std::string, std::string> &potential,
                                Real x, Real R, Real delta, DomainSide side)
{
    return case3WindowFunction(x, R, delta, side) * evaluateFunction(potential, x);
}

AsymptoteClassification classifyAsymptote(const std::map<std::string, std::string> &potential,
                                           const SpatialDomain &domain,
                                           DomainSide side,
                                           std::ostream &warnOut)
{
    const Real reference = (side == DomainSide::Left) ? domain.xMin : domain.xMax;
    const Real sign = (side == DomainSide::Left) ? -1.0 : 1.0;
    const Real scale = std::max(std::abs(reference), 1.0);

    constexpr int kNumSamples = 16;
    constexpr Real kRatio = 4.0;
    std::vector<Real> V(kNumSamples);
    for (int k = 0; k < kNumSamples; ++k)
        V[k] = evaluateFunction(potential, reference + sign * scale * std::pow(kRatio, k));

    const ConvergenceFit fit = classifySequenceConvergence(V, kRatio);
    const Real nan = std::numeric_limits<Real>::quiet_NaN();

    AsymptoteClassification result{};
    result.warningEmitted = false;

    if (fit.isDivergent)
    {
        result.asymptoteCase = AsymptoteCase::HardWall;
        result.subType = AsymptoteSubType::NotApplicable;
        result.fittedAsymptoticValue = nan;
        result.powerLawExponent = nan;
        result.recommendedTransitionWidth = nan;
        return result;
    }

    if (fit.isFlat)
    {
        result.asymptoteCase = AsymptoteCase::AnalyticAsymptote;
        result.subType = AsymptoteSubType::Flat;
        result.fittedAsymptoticValue = fit.fittedLimit;
        result.powerLawExponent = nan;
        result.recommendedTransitionWidth = nan;
        return result;
    }

    if (std::abs(fit.powerLawExponent - 1.0) <= 0.15)
    {
        result.asymptoteCase = AsymptoteCase::AnalyticAsymptote;
        result.subType = AsymptoteSubType::Coulomb;
        result.fittedAsymptoticValue = fit.fittedLimit;
        result.powerLawExponent = fit.powerLawExponent;
        result.recommendedTransitionWidth = nan;
        return result;
    }

    result.asymptoteCase = AsymptoteCase::Irregular;
    result.subType = AsymptoteSubType::NotApplicable;
    result.fittedAsymptoticValue = fit.fittedLimit;
    result.powerLawExponent = fit.powerLawExponent;
    result.recommendedTransitionWidth = 0.1 * (domain.xMax - domain.xMin);
    result.warningEmitted = true;

    warnOut << "Warning: potential asymptote on the " << (side == DomainSide::Left ? "left" : "right")
            << " side is irregular (fitted power-law exponent p=" << fit.powerLawExponent
            << "); the potential will be smoothly tapered to zero over a transition width delta="
            << result.recommendedTransitionWidth
            << " approaching the box boundary (see docs/planning/boundary-condition-case-3-smoothing.md), "
            << "avoiding an abrupt truncation. This remains an approximation -- the true asymptotic tail "
            << "is not analytically known -- so continuum normalization will be approximate.\n";

    return result;
}

std::pair<std::vector<Real>, std::vector<Real>>
fillBandedMatrices(const bspline::BSpline &bs, int nEn, int order, int L, std::map<std::string, std::string> potential)
{
    std::vector<Real> Hmat(order * nEn, 0.0);
    std::vector<Real> Smat(order * nEn, 0.0);

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

EigenResult solveTISE(int nNodes, int order, Real rMin, Real rMax, int L, std::map<std::string, std::string> potential)
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
    return solveGeneralizedEigenproblem(std::move(H), std::move(S), nEn, order);
}

} // namespace tise
