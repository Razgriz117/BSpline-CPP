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
//
// Regex/bound-parsing logic factored into parseInterval() (below) so A4a's
// detectPotentialStructure can reuse it to get numeric piece bounds without
// duplicating this nontrivial regex. File-local only -- not exposed in
// tise.hpp, since no caller outside this file needs the numeric bounds
// directly.
namespace
{
struct ParsedInterval
{
    double lower;
    double upper;
    bool lowerInclusive;
    bool upperInclusive;
};

ParsedInterval parseInterval(const std::string& interval)
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

    ParsedInterval out;
    out.lowerInclusive = m[1] == "[";
    out.upperInclusive = m[4] == "]";
    out.lower = parseBound(m[2]);
    out.upper = parseBound(m[3]);
    return out;
}
} // namespace

bool inInterval(double x, const std::string& interval)
{
    const ParsedInterval iv = parseInterval(interval);
    bool leftOK = iv.lowerInclusive ? (x >= iv.lower) : (x > iv.lower);
    bool rightOK = iv.upperInclusive ? (x <= iv.upper) : (x < iv.upper);
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

BoundStateClassification classifyBoundStates(const EigenResult &result, Real threshold)
{
    BoundStateClassification out;
    out.isBound.resize(result.values.size());
    out.nBound = 0;
    for (std::size_t i = 0; i < result.values.size(); ++i)
    {
        bool bound = result.values[i] < threshold;
        out.isBound[i] = bound;
        if (bound)
            ++out.nBound;
    }
    return out;
}

// === A3: well-containment diagnostic (SDD Sec. 5.2.3, Sec. 6.4) ===
ContainmentCheck checkWellContainment(const bspline::BSpline &bs,
                                       const std::vector<Real> &coeffs,
                                       Real xBoundary,
                                       Real tol)
{
    const int n = static_cast<int>(coeffs.size());
    ContainmentCheck out;
    out.psiPrimeAtBoundary = bs.eval(xBoundary, coeffs.data(), n, 1);
    out.notWellContained = std::abs(out.psiPrimeAtBoundary) > tol;
    return out;
}

// === A4a: strategic node placement (REQ-F-050), required core only ===
namespace
{
// Shrinking-offset probe toward x0 from one side, reusing A1's numeric core
// (classifySequenceConvergence) but inverted: offsets shrink toward x0
// instead of growing toward infinity -- valid because that function's
// divergence test only compares finite-difference magnitudes across the
// given array, with no built-in assumption about sample direction (see
// docs/planning/engineer-a-plan-A4.md). `sign` is -1 to sample the piece
// left of x0, +1 for the piece right of x0; `pieceWidth` is that piece's
// own width, used to clamp the starting offset so the probe can never step
// outside a narrow piece (e.g. the box+barrier's width-1 barrier slab) --
// A1's own scale formula alone isn't safe for this inward direction.
bool isSingularApproaching(const std::map<std::string, std::string> &potential,
                            Real x0, Real sign, Real pieceWidth)
{
    constexpr int kNumSamples = 16;
    constexpr Real kRatio = 4.0;
    const Real scale = std::max(std::abs(x0), 1.0);
    const Real h0 = std::isfinite(pieceWidth) ? std::min(scale, 0.25 * pieceWidth) : scale;

    std::vector<Real> V(kNumSamples);
    for (int k = 0; k < kNumSamples; ++k)
        V[k] = evaluateFunction(potential, x0 + sign * h0 / std::pow(kRatio, k));

    return classifySequenceConvergence(V, kRatio).isDivergent;
}

// One-sided value estimate at x0, approaching from `sign` (-1 = left piece,
// +1 = right piece). Richardson-extrapolated over step sizes h and h/2 to
// cancel the O(h) bias from necessarily sampling strictly inside the piece
// (open piece boundaries, e.g. "[0,5)", never include x0=5 itself):
// V(x0+sign*h) = V0 + sign*h*V0' + O(h^2), so 2*V(h/2) - V(h) = V0 + O(h^2).
Real oneSidedValue(const std::map<std::string, std::string> &potential,
                    Real x0, Real sign, Real h)
{
    const Real gH  = evaluateFunction(potential, x0 + sign * h);
    const Real gH2 = evaluateFunction(potential, x0 + sign * h / 2.0);
    return 2.0 * gH2 - gH;
}

// One-sided slope (dV/dx) estimate at x0, same side convention. Builds a
// standard one-sided 2nd-order finite-difference derivative,
// f'(a) ~ [-3f(a) + 4f(a+h) - f(a+2h)] / (2h), centered at a = x0+sign*h
// (not x0 itself, since x0 may not be in this piece's domain) -- that
// off-centering is an O(h) bias, which the same Richardson trick (evaluate
// at h and h/2, then 2*D(h/2) - D(h)) cancels, leaving O(h^2) error --
// negligible at h~1e-4 against the StitchedKink threshold below.
Real oneSidedSlope(const std::map<std::string, std::string> &potential,
                    Real x0, Real sign, Real h)
{
    auto forwardDeriv = [&](Real step) {
        const Real g1 = evaluateFunction(potential, x0 + sign * step);
        const Real g2 = evaluateFunction(potential, x0 + sign * 2.0 * step);
        const Real g3 = evaluateFunction(potential, x0 + sign * 3.0 * step);
        return sign * (-3.0 * g1 + 4.0 * g2 - g3) / (2.0 * step);
    };
    const Real dH  = forwardDeriv(h);
    const Real dH2 = forwardDeriv(h / 2.0);
    return 2.0 * dH2 - dH;
}
} // namespace

std::vector<DetectedJoin> detectPotentialStructure(const std::map<std::string, std::string> &potential)
{
    struct PieceInfo { Real lower, upper; };
    std::vector<PieceInfo> pieces;
    pieces.reserve(potential.size());
    for (const auto &entry : potential)
    {
        const ParsedInterval iv = parseInterval(entry.first);
        pieces.push_back({iv.lower, iv.upper});
    }
    // Sort by numeric lower bound -- std::map's native order is lexicographic
    // on the domain *string* (e.g. "[10,...)" sorts before "[5,...)"), which
    // would silently misorder joins if used directly.
    std::sort(pieces.begin(), pieces.end(),
              [](const PieceInfo &a, const PieceInfo &b) { return a.lower < b.lower; });

    constexpr Real kValueJumpTol = 1e-6; // O(1)-O(10) potential magnitudes in
                                          // this codebase; a genuine Step is
                                          // an O(1) jump, ~6 orders of margin.
    constexpr Real kSlopeJumpTol = 1e-4; // one-sided slope estimate's O(h^2)
                                          // truncation error is ~1e-8 at
                                          // h=1e-4, ~4 orders below this.
    constexpr Real kFdStep = 1e-4;

    auto clampedStep = [&](Real width) {
        return std::isfinite(width) ? std::min(kFdStep, 0.1 * width) : kFdStep;
    };

    auto classifyJoin = [&](Real x0, Real leftWidth, Real rightWidth) {
        if (isSingularApproaching(potential, x0, -1.0, leftWidth) ||
            isSingularApproaching(potential, x0, +1.0, rightWidth))
            return JoinType::Singular;

        const Real hL = clampedStep(leftWidth);
        const Real hR = clampedStep(rightWidth);
        const Real vL = oneSidedValue(potential, x0, -1.0, hL);
        const Real vR = oneSidedValue(potential, x0, +1.0, hR);
        if (std::abs(vL - vR) > kValueJumpTol)
            return JoinType::Step;

        const Real sL = oneSidedSlope(potential, x0, -1.0, hL);
        const Real sR = oneSidedSlope(potential, x0, +1.0, hR);
        if (std::abs(sL - sR) > kSlopeJumpTol)
            return JoinType::StitchedKink;

        return JoinType::Continuous;
    };

    std::vector<DetectedJoin> joins;

    for (std::size_t i = 0; i + 1 < pieces.size(); ++i)
    {
        const Real x0 = pieces[i].upper;
        const Real leftWidth = pieces[i].upper - pieces[i].lower;
        const Real rightWidth = pieces[i + 1].upper - pieces[i + 1].lower;
        joins.push_back({x0, classifyJoin(x0, leftWidth, rightWidth)});
    }

    // Global domain edges: only flagged if singular (an ordinary box wall
    // isn't interesting to strategic node placement). Also probes each
    // piece's own endpoints, not just inter-piece joins, so a single-piece
    // potential like {"(0, inf)": "-1/x + 1/x^2"} still gets its x=0
    // singularity flagged even though it has no "join" at all.
    if (!pieces.empty())
    {
        const PieceInfo &first = pieces.front();
        if (std::isfinite(first.lower))
        {
            const Real width = first.upper - first.lower;
            if (isSingularApproaching(potential, first.lower, +1.0, width))
                joins.push_back({first.lower, JoinType::Singular});
        }
        const PieceInfo &last = pieces.back();
        if (std::isfinite(last.upper))
        {
            const Real width = last.upper - last.lower;
            if (isSingularApproaching(potential, last.upper, -1.0, width))
                joins.push_back({last.upper, JoinType::Singular});
        }
    }

    return joins;
}

std::vector<StrategicKnot> strategicKnotsFromJoins(const std::vector<DetectedJoin> &joins, int order)
{
    std::vector<StrategicKnot> knots;
    knots.reserve(joins.size());
    for (const auto &j : joins)
    {
        int extra = 0;
        switch (j.type)
        {
            case JoinType::Step:         extra = std::max(0, order - 3); break;
            case JoinType::StitchedKink: extra = std::max(0, order - 4); break;
            case JoinType::Singular:
            case JoinType::Continuous:   extra = 0; break;
        }
        if (extra > 0)
            knots.push_back({j.x, extra});
    }
    return knots;
}

std::vector<Real> buildStrategicRadialGrid(int nNodes, Real rMin, Real rMax,
                                            const std::vector<StrategicKnot> &knots)
{
    std::vector<Real> grid = buildUniformRadialGrid(nNodes, rMin, rMax);

    std::vector<StrategicKnot> sorted = knots;
    std::sort(sorted.begin(), sorted.end(),
              [](const StrategicKnot &a, const StrategicKnot &b) { return a.x < b.x; });

    for (const auto &knot : sorted)
    {
        auto it = std::lower_bound(grid.begin(), grid.end(), knot.x);
        const bool alreadyPresent = (it != grid.end()) && (std::abs(*it - knot.x) < 1e-12);

        if (!alreadyPresent)
            it = grid.insert(it, knot.x); // splice a new base point at knot.x

        grid.insert(it + 1, knot.extraMultiplicity, knot.x); // extra degenerate copies
    }

    return grid;
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
