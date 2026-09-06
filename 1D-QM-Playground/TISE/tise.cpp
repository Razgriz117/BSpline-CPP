#include "tise.hpp"

#define _USE_MATH_DEFINES

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
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

void validateNoOverlappingPotentialPieces(const std::map<std::string, std::string> &potential)
{
    std::vector<std::pair<std::string, ParsedInterval>> pieces;
    pieces.reserve(potential.size());
    for (const auto &[domain, fn] : potential)
    {
        (void)fn;
        pieces.emplace_back(domain, parseInterval(domain));
    }

    // A large-but-finite stand-in for +/-infinity, usable as a concrete
    // evaluateFunction-style probe x -- mirrors the same "probe far beyond
    // an infinite bound" idea tise_solver_main.cpp's own Case-3 detection
    // uses, though with its own separate constant (each solves a distinct
    // problem -- this one just needs *some* finite number outside every
    // finite piece to test interval membership against; unifying the two
    // idioms is not worth doing for that alone).
    // kUnboundedProbeMagnitude divides by 4 (rather than using max() itself)
    // to leave headroom against overflow in later arithmetic on the probe
    // value (e.g. this function's own midpoint averaging below, or a
    // caller doing further arithmetic on the reported x in its error
    // message) -- max() itself has no room left to be added to or doubled
    // without overflowing to infinity.
    constexpr double kUnboundedProbeMagnitude = std::numeric_limits<double>::max() / 4.0;
    auto finiteProbe = [](double bound, bool isUpper) {
        if (std::isinf(bound))
            return isUpper ? kUnboundedProbeMagnitude : -kUnboundedProbeMagnitude;
        return bound;
    };

    for (std::size_t i = 0; i < pieces.size(); ++i)
    {
        for (std::size_t j = i + 1; j < pieces.size(); ++j)
        {
            const auto &[domainA, ivA] = pieces[i];
            const auto &[domainB, ivB] = pieces[j];

            std::vector<double> probes = {
                finiteProbe(ivA.lower, false), finiteProbe(ivA.upper, true),
                finiteProbe(ivB.lower, false), finiteProbe(ivB.upper, true),
            };
            if (std::isfinite(ivA.lower) && std::isfinite(ivA.upper))
                probes.push_back((ivA.lower + ivA.upper) / 2.0);
            if (std::isfinite(ivB.lower) && std::isfinite(ivB.upper))
                probes.push_back((ivB.lower + ivB.upper) / 2.0);

            for (double x : probes)
            {
                if (inInterval(x, domainA) && inInterval(x, domainB))
                    throw std::runtime_error("overlapping potential piece domains: '" + domainA +
                                              "' and '" + domainB + "' both cover x=" + std::to_string(x));
            }
        }
    }
}

std::pair<std::string, std::string> parsePotentialPiece(const std::string &piece)
{
    auto first = piece.find_first_not_of(" \t");
    auto last  = piece.find_last_not_of(" \t");
    if (first == std::string::npos || piece[first] != '{' || piece[last] != '}')
        throw std::runtime_error("malformed potential piece (expected a {...} dict literal): " + piece);

    static const std::regex domainRe(R"('domain'\s*:\s*'([^']*)')");
    static const std::regex functionRe(R"('function'\s*:\s*'([^']*)')");

    std::smatch mDomain, mFunction;
    if (!std::regex_search(piece, mDomain, domainRe))
        throw std::runtime_error("potential piece missing 'domain' key: " + piece);
    if (!std::regex_search(piece, mFunction, functionRe))
        throw std::runtime_error("potential piece missing 'function' key: " + piece);

    return {mDomain[1].str(), mFunction[1].str()};
}

// === A1: boundary-condition asymptote classifier (REQ-F-030) ===
// See classifySequenceConvergence's own declaration comment in tise.hpp for
// what this function decides; the thresholds below are explained where each
// is used.
ConvergenceFit classifySequenceConvergence(const std::vector<Real> &V, Real ratio,
                                            Real flatnessAbsTol, Real flatnessRelTol,
                                            Real shrinkFactorFloor, Real divergenceThreshold,
                                            int windowSize)
{
    // The tail-window power-law fit below (reached when V is neither flat
    // nor divergent) needs at least one successive-difference-ratio
    // estimate: windowStart = max(0, N-windowSize) and windowEnd = N-3 must
    // satisfy windowStart <= windowEnd, which reduces to exactly N >= 3 AND
    // windowSize >= 3 (checked independently: e.g. N=100 with windowSize=2
    // still gives windowStart=N-2 > windowEnd=N-3, no valid k; N=2 with
    // windowSize=100 still gives windowEnd=-1, no valid k). Below that
    // minimum, pEstimates would end up empty and pEstimates[size()/2] --
    // this function's own median step -- would read out of bounds. This was
    // unreachable while windowSize was a hardcoded literal (7 >= 3); it
    // became caller-reachable once windowSize turned into a parameter, so
    // it's validated explicitly now rather than left as a latent UB trap.
    if (V.size() < 3)
        throw std::runtime_error("classifySequenceConvergence: V must have at least 3 samples");
    if (windowSize < 3)
        throw std::runtime_error("classifySequenceConvergence: windowSize must be >= 3");

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

    // Flatness pre-check: successive differences are all numerically zero,
    // to within an absolute tolerance (flatnessAbsTol) plus a tolerance
    // relative to the sequence's own magnitude (flatnessRelTol * maxAbsV) --
    // the relative term keeps this check meaningful for a V that is large in
    // magnitude but still flat, not just one near zero.
    if (maxAbsDV <= flatnessAbsTol + flatnessRelTol * maxAbsV)
    {
        fit.isDivergent = false;
        fit.isFlat = true;
        fit.fittedLimit = V[N - 1];
        fit.powerLawExponent = nan;
        return fit;
    }

    // Divergence check: compare difference magnitude at the very end of the
    // window vs. the middle (shrinkFactorFloor guards the denominator
    // against division by a near-zero difference). If the end-to-middle
    // ratio hasn't shrunk below divergenceThreshold, successive differences
    // aren't shrinking fast enough for the sequence to be settling toward a
    // finite limit -- treat it as divergent (e.g. approaching a 1/x pole).
    const int mid = (N - 2) / 2;
    const Real shrinkFactor = std::abs(dV[N - 2]) / std::max(std::abs(dV[mid]), shrinkFactorFloor);
    if (shrinkFactor >= divergenceThreshold)
    {
        fit.isDivergent = true;
        fit.isFlat = false;
        fit.fittedLimit = nan;
        fit.powerLawExponent = nan;
        return fit;
    }

    // Power-law fit via successive-ratio over a tail window, avoiding both
    // near-field transients (early samples, still settling) and far-field
    // floating-point noise (very late samples, where dV is tiny and
    // dominated by rounding error). The window starts `windowSize` samples
    // back from the end of V (windowStart = N - windowSize); the "-3" below
    // is unrelated to windowSize -- it's fixed headroom so the k+1 lookahead
    // a few lines down never reads past the last valid finite difference,
    // regardless of windowSize. Net effect: the fit uses windowSize-2
    // successive-difference-ratio estimates (e.g. windowSize=7 -> 5
    // estimates), not windowSize differences directly.
    const int windowStart = std::max(0, N - windowSize);
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
// kDefaultAsymptoteNumSamples/kDefaultAsymptoteRatio/kDefaultAsymptoteBackupScale
// (shared by classifyAsymptote below and isSingularApproaching further down
// this file) live in tise.hpp, not here -- see the comment there for why:
// classifyAsymptote's own default-argument list in the header must be able
// to name them, which an anonymous-namespace constant confined to this .cpp
// file cannot provide.

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
                                           std::ostream &warnOut,
                                           int numSamples,
                                           Real ratio,
                                           Real backupScale,
                                           Real coulombExponentTol,
                                           Real transitionWidthFraction)
{
    // numSamples feeds classifySequenceConvergence's V directly (below), and
    // that function's own tail-window power-law fit needs at least 3 samples
    // to produce a single successive-difference-ratio estimate (see its own
    // guard for the derivation) -- fail fast here, before doing any sampling
    // work or constructing V, so a bad numSamples reports against this
    // function's own name rather than surfacing as a confusing failure deep
    // inside classifySequenceConvergence.
    if (numSamples < 3)
        throw std::runtime_error("classifyAsymptote: numSamples must be >= 3");

    const Real reference = (side == DomainSide::Left) ? domain.xMin : domain.xMax;
    const Real sign = (side == DomainSide::Left) ? -1.0 : 1.0;
    // scale is the offset of the first sample from `reference`; backupScale
    // floors it away from zero when the domain boundary itself is at/near
    // x=0 (e.g. the left edge of a radial [0, rMax) domain), where
    // |reference| alone would otherwise produce a degenerate zero-width
    // first step.
    const Real scale = std::max(std::abs(reference), backupScale);

    std::vector<Real> V(numSamples);
    for (int k = 0; k < numSamples; ++k)
        V[k] = evaluateFunction(potential, reference + sign * scale * std::pow(ratio, k));

    const ConvergenceFit fit = classifySequenceConvergence(V, ratio);
    const Real nan = std::numeric_limits<Real>::quiet_NaN();

    AsymptoteClassification result{};
    result.warningEmitted = false;

    if (fit.isDivergent)
    {
        result.asymptoteCase = AsymptoteCase::HardWall;
        result.subType = AsymptoteSubType::NotApplicable;
        result.fittedAsymptoticValue = nan;
        result.powerLawExponent = nan;
        result.fittedPowerLawCoefficient = nan;
        result.recommendedTransitionWidth = nan;
        return result;
    }

    if (fit.isFlat)
    {
        result.asymptoteCase = AsymptoteCase::AnalyticAsymptote;
        result.subType = AsymptoteSubType::Flat;
        result.fittedAsymptoticValue = fit.fittedLimit;
        result.powerLawExponent = nan;
        result.fittedPowerLawCoefficient = nan;
        result.recommendedTransitionWidth = nan;
        return result;
    }

    // Coefficient C in V ~ fittedLimit + C/x^p, read off the last (farthest,
    // least noisy relative to the fitted tail) sample: V[last] = fittedLimit
    // + C/|x_last|^p => C = (V[last]-fittedLimit)*|x_last|^p. abs() on
    // x_last: the physical tail falloff is naturally expressed in distance
    // from the origin, not signed coordinate (matters for DomainSide::Left,
    // where samples run negative).
    const Real xLast = reference + sign * scale * std::pow(ratio, numSamples - 1);
    const Real fittedCoefficient = (V[numSamples - 1] - fit.fittedLimit) *
                                    std::pow(std::abs(xLast), fit.powerLawExponent);

    // A true Coulomb tail (V ~ C/r) is a power law with exponent exactly 1;
    // coulombExponentTol is how far the fitted exponent may drift from 1
    // (numerical fitting noise, finite sampling) and still be called
    // Coulomb rather than an unrecognized/irregular power law.
    if (std::abs(fit.powerLawExponent - 1.0) <= coulombExponentTol)
    {
        result.asymptoteCase = AsymptoteCase::AnalyticAsymptote;
        result.subType = AsymptoteSubType::Coulomb;
        result.fittedAsymptoticValue = fit.fittedLimit;
        result.powerLawExponent = fit.powerLawExponent;
        result.fittedPowerLawCoefficient = fittedCoefficient;
        result.recommendedTransitionWidth = nan;
        return result;
    }

    result.asymptoteCase = AsymptoteCase::Irregular;
    result.subType = AsymptoteSubType::NotApplicable;
    result.fittedAsymptoticValue = fit.fittedLimit;
    result.powerLawExponent = fit.powerLawExponent;
    result.fittedPowerLawCoefficient = fittedCoefficient;
    // Neither flat nor Coulomb -- the tail shape is unrecognized. Recommend
    // tapering it to zero (via case3WindowFunction) over a transition width
    // that scales with the box size (transitionWidthFraction of the total
    // domain span), so the taper is neither vanishingly thin (would still
    // look like an abrupt truncation to the basis) nor so wide it eats into
    // the trusted interior.
    result.recommendedTransitionWidth = transitionWidthFraction * (domain.xMax - domain.xMin);
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

namespace
{
// colOf[i] = 1-based matrix column for physical B-spline i (1-based), or 0
// if i is excluded by dropSet. Size nBSplines+1 (index 0 unused). nKept is
// the true number of surviving indices.
struct ColumnMap
{
    std::vector<int> colOf;      // colOf[physicalIdx] = column (1-based), 0 if dropped
    std::vector<int> physicalOf; // physicalOf[col] = physicalIdx (1-based); size nBSplines+1,
                                  // only entries 1..nKept meaningful
    int nKept;
};

ColumnMap columnIndexMap(int nBSplines, const std::vector<int> &dropSet)
{
    std::vector<bool> isDropped(nBSplines + 1, false);
    for (int d : dropSet)
    {
        if (d < 1 || d > nBSplines)
            throw std::runtime_error("dropSet index out of range [1," +
                                      std::to_string(nBSplines) + "]: " + std::to_string(d));
        isDropped[d] = true;
    }
    std::vector<int> colOf(nBSplines + 1, 0);
    std::vector<int> physicalOf(nBSplines + 1, 0);
    int c = 0;
    for (int idx = 1; idx <= nBSplines; ++idx)
        if (!isDropped[idx])
        {
            colOf[idx] = ++c;
            physicalOf[c] = idx;
        }
    return {colOf, physicalOf, c};
}

// Shared by fillBandedMatrices/eigenstateCoefficients: resolves the
// caller's optional drop-set to the classic {1, nBSplines} default when
// absent, and validates nEn matches the true kept count. Both functions
// must reach this exact same check -- a mismatch here previously risked a
// silent out-of-bounds write: e.g. nBSplines=10, order=4, dropSet={1,10}
// leaves 8 kept indices, but if a caller separately passed nEn=7, the band
// matrices below are allocated order*nEn = 4*7 elements while the column
// map (built from dropSet alone, independent of nEn) still produces column
// indices up to 8 -- the last kept column's write lands a full `order`
// (=4) elements past the end of the buffer, corrupting adjacent memory
// with no bounds check to catch it. Checking nKept == nEn here, once, up
// front, turns that into an immediate std::runtime_error instead.
ColumnMap resolveDropSet(int nBSplines, int nEn, const std::optional<std::vector<int>> &dropSet)
{
    std::vector<int> drop = dropSet.value_or(std::vector<int>{1, nBSplines});
    ColumnMap map = columnIndexMap(nBSplines, drop);
    if (map.nKept != nEn)
        throw std::runtime_error("nEn (" + std::to_string(nEn) +
                                  ") does not match nBSplines - |dropSet| (" +
                                  std::to_string(map.nKept) + ")");
    return map;
}

// Shared by matchAsymptotic/writeContinuumInfo: transform one energy's
// buildContinuumState coefficients (in the confined-eigenstate basis
// {phi_n} plus B_N) into true B-spline coefficients, suitable for
// bs.eval. `stateCoeffs` is one buildContinuumState row (size nEn+1,
// last entry always 1.0, the B_N coefficient); `physicalOf` is
// columnIndexMap's output for whatever drop-set produced `eigen`.
// Deliberately factored out after this exact transform was found
// duplicated-and-diverged once already (writeContinuumInfo had the
// pre-A4b, un-transformed version while matchAsymptotic had the fix) --
// keeping one copy means a future correction can't reintroduce that bug.
std::vector<Real> continuumStateToBSplineCoeffs(const std::vector<Real> &stateCoeffs,
                                                 const EigenResult &eigen,
                                                 int nBSplines,
                                                 const std::vector<int> &physicalOf)
{
    const int nEn = eigen.dim;
    std::vector<Real> fc(nBSplines, 0.0);
    for (int j = 0; j < nEn; ++j)
    {
        Real coeff = 0.0;
        for (int n = 0; n < nEn; ++n)
            coeff += stateCoeffs[n] * eigen.vectors[n * eigen.ldz + j];
        fc[physicalOf[j + 1] - 1] = coeff;
    }
    fc[nBSplines - 1] = stateCoeffs[nEn];
    return fc;
}
} // namespace

std::pair<std::vector<Real>, std::vector<Real>>
fillBandedMatrices(const bspline::BSpline &bs, int nEn, int order, int L,
                    std::map<std::string, std::string> potential,
                    std::optional<std::vector<int>> dropSet,
                    std::optional<Real> case3RightR,
                    std::optional<Real> case3RightDelta)
{
    const int nBSplines = bs.getNBSplines();
    const ColumnMap map = resolveDropSet(nBSplines, nEn, dropSet);
    const std::vector<int> &colOf = map.colOf;

    std::vector<Real> Hmat(order * nEn, 0.0);
    std::vector<Real> Smat(order * nEn, 0.0);

    bspline::D2DFun fUni = [](double, const double *) { return 1.0; };
    // Piecewise potential supplied by the caller, evaluated per-x via muparser.
    // `L` is no longer used to select the potential here; it is retained for
    // eigenvalueError()'s comparison against the analytic hydrogen spectrum.
    // Case-3 remediation: when the caller supplies both case3RightR/Delta
    // (from classifyAsymptote's Irregular classification), evaluate through
    // the windowed/tapered potential instead of the raw one.
    bspline::D2DFun fPot;
    if (case3RightR.has_value() && case3RightDelta.has_value())
    {
        const Real R = *case3RightR;
        const Real delta = *case3RightDelta;
        fPot = [potential, R, delta](double x, const double *) {
            return evaluateWindowedPotential(potential, x, R, delta, DomainSide::Right);
        };
    }
    else
    {
        fPot = [potential](double x, const double *) {
            return evaluateFunction(potential, x);
        };
    }
    // Previous hardcoded radial hydrogen-like potential, kept for reference:
    // bspline::D2DFun fPot = [L](double x, const double *) {
    //     return radialPotential(x, L);
    // };
    double parvec[1] = {0.0};

    auto bandIndex = [&](int row, int col) {
        return (row - 1) + (col - 1) * order;
    };

    for (int iBs2 = 1; iBs2 <= nBSplines; ++iBs2)
    {
        const int col2 = colOf[iBs2];
        if (col2 == 0)
            continue;
        const int iBs1Min = std::max(1, iBs2 - order + 1);
        for (int iBs1 = iBs1Min; iBs1 <= iBs2; ++iBs1)
        {
            const int col1 = colOf[iBs1];
            if (col1 == 0)
                continue;

            Real overlap       = bs.integral(fUni, iBs1, iBs2);
            Real kinetic       = bs.integral(fUni, iBs1, iBs2, 1, 1) / 2.0;
            Real potentialTerm = bs.integral(fPot, iBs1, iBs2, 0, 0, parvec);

            const int row = order + col1 - col2;
            const int idx = bandIndex(row, col2);

            Smat[idx] = overlap;
            Hmat[idx] = kinetic + potentialTerm;
        }
    }

    return {Hmat, Smat};
}

std::pair<std::vector<Real>, std::vector<Real>> precomputeBoundaryCoupling(
    int order, int nEn,
    std::vector<Real> Hmat,
    std::vector<Real> Smat,
    EigenResult eigen,
    std::optional<int> nBSplinesOpt,
    std::optional<std::vector<int>> dropSet)
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
    //
    // Generalized via colOf/physicalOf (A4b, REQ-F-050) instead of the
    // classic hardcoded "only B_1 dropped" shift-by-1 arithmetic: when
    // exactly B_1 is dropped, physical index i always maps to matrix
    // column i-1, so the old code could get away with a plain "-1" offset
    // everywhere. Once solveTISE can also drop interior B-splines (A4b
    // singular-join removal), that one-to-one offset no longer holds --
    // physical index i's column depends on how many OTHER indices below it
    // were also dropped -- so every physical-to-column lookup here goes
    // through colOf/physicalOf instead of assuming a fixed shift.
    const int nBSplines = nBSplinesOpt.value_or(nEn + 2);
    const ColumnMap map = columnIndexMap(nBSplines, dropSet.value_or(std::vector<int>{1}));
    const std::vector<int> &colOf = map.colOf;
    const int iBs2 = nBSplines; // B_N, always kept (never in dropSet)
    const int col2 = colOf[iBs2];
    if (map.nKept != nEn + 1 || col2 != nEn + 1)
        throw std::runtime_error("precomputeBoundaryCoupling: nEn inconsistent with nBSplines/dropSet "
                                  "(expected nBSplines - |dropSet| == nEn + 1, and B_N must be kept)");

    std::vector<Real> HB_N(nEn, 0.0), SB_N(nEn, 0.0);
    const int iBs1Min = std::max(1, iBs2 - order + 1);
    for (int iBs1 = iBs1Min; iBs1 <= iBs2 - 1; ++iBs1)
    {
        const int col1 = colOf[iBs1];
        if (col1 == 0)
            continue;
        const int row = order + col1 - col2;      // matches fillBandedMatrices' bandIndex row formula
        const int idx = (row - 1) + (col2 - 1) * order;
        HB_N[col1 - 1] = Hmat[idx];
        SB_N[col1 - 1] = Smat[idx];
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
    std::vector<Real> grid,
    std::optional<int> nBSplinesOpt,
    std::optional<std::vector<int>> dropSet,
    Real poleTolFraction,
    std::ostream &warnOut)
{
    auto [coeffs1, coeffs2] = precomputeBoundaryCoupling(order, nEn, Hmat, Smat, eigen, nBSplinesOpt, dropSet);

    std::vector<std::vector<Real>> states(grid.size(), std::vector<Real>(eigen.values.size() + 1, 0.0));

    const Real poleTol = grid.size() >= 2 ? poleTolFraction * (grid[1] - grid[0]) : Real(0);

    // Compute the coeffs for each point on the energy grid
    for(int E_idx = 0; E_idx < grid.size(); ++E_idx)
    {
        for(int i = 0; i < nEn; ++i)
        {
            states[E_idx][i] = (coeffs1[i] - (grid[E_idx] * coeffs2[i])) / (grid[E_idx] - eigen.values[i]);
        }
        states[E_idx][nEn] = 1.0;

        if (poleTol > 0.0 && nEn > 0)
        {
            int closest = 0;
            Real closestGap = std::abs(grid[E_idx] - eigen.values[0]);
            for (int i = 1; i < nEn; ++i)
            {
                Real gap = std::abs(grid[E_idx] - eigen.values[i]);
                if (gap < closestGap)
                {
                    closestGap = gap;
                    closest = i;
                }
            }
            if (closestGap < poleTol)
                warnOut << "Warning: continuum energy grid point E=" << grid[E_idx]
                        << " is within " << closestGap << " of confined eigenvalue E_" << closest
                        << "=" << eigen.values[closest] << "; this energy's continuum state is "
                        << "likely a finite-box discretization artifact (see ADR-0007), not a "
                        << "physical feature -- treat with suspicion.\n";
        }
    }

    return states;

}

// === Coulomb-tail continuum matching (ADR-0009, supersedes ADR-0010) ===
// See the module-level comment on matchAsymptotic's declaration (tise.hpp)
// and docs/planning/coulomb-tail-continuum-matching.md for the full
// derivation and validation record.

namespace
{

// Lanczos approximation for the complex Gamma function, g=7/n=9 coefficient
// set (a standard, widely-published choice -- e.g. the one used by Numerical
// Recipes' own complex-Gamma routine). Only the ARGUMENT of Gamma(1+i*eta)
// is actually needed (coulombPhaseShift below), but implementing the full
// complex Gamma is no harder than a real-argument-only special case and is
// easier to validate directly against a reference (mpmath.gamma). Validated
// during implementation: arg(lanczosGamma(1+i*eta)) matches mpmath's
// arg(gamma(1+i*eta)) to ~1e-15 across eta in [-5, 5].
std::complex<Real> lanczosGamma(std::complex<Real> z)
{
    static constexpr int kG = 7;
    static constexpr Real kCoeffs[kG + 2] = {
        0.99999999999980993,   676.5203681218851,     -1259.1392167224028,
        771.32342877765313,    -176.61502916214059,   12.507343278686905,
        -0.13857109526572012,  9.9843695780195716e-6, 1.5056327351493116e-7,
    };

    if (z.real() < 0.5)
        // Reflection formula: extends the approximation (valid for
        // Re(z)>=0.5) to the rest of the complex plane. Only Re(z)=1
        // (from coulombPhaseShift's Gamma(1+i*eta)) is actually exercised
        // today, so this branch is untested dead code in practice, but
        // it's cheap and standard to include for a self-contained,
        // correct-for-any-z implementation.
        return kPi / (std::sin(kPi * z) * lanczosGamma(Real(1.0) - z));

    z -= Real(1.0);
    std::complex<Real> x(kCoeffs[0], 0.0);
    for (int i = 1; i < kG + 2; ++i)
        x += kCoeffs[i] / (z + Real(i));
    std::complex<Real> t = z + Real(kG) + Real(0.5);
    return std::sqrt(2.0 * kPi) * std::pow(t, z + Real(0.5)) * std::exp(-t) * x;
}

} // namespace

Real coulombPhaseShift(int l, Real eta)
{
    Real sigma = std::arg(lanczosGamma(std::complex<Real>(1.0, eta)));
    for (int k = 1; k <= l; ++k)
        sigma += std::atan(eta / k);
    return sigma;
}

CoulombWaveResult evaluateCoulombFunctions(int l, Real eta, Real rho,
                                            Real farMultiplier, Real stepSize)
{
    // Local "index of refraction" for the Coulomb radial equation
    // u'' + k(rho)^2 u = 0, k(rho)^2 = 1 - 2*eta/rho - l(l+1)/rho^2.
    auto localK = [&](Real r) {
        return std::sqrt(1.0 - 2.0 * eta / r - l * (l + 1) / (r * r));
    };
    // Central-difference derivative of localK -- cheap, and only ever
    // evaluated once per F/G pair (at rhoStart), not per integration step.
    auto localKPrime = [&](Real r) {
        constexpr Real h = 1e-4;
        return (localK(r + h) - localK(r - h)) / (2.0 * h);
    };

    const Real sigmaL = coulombPhaseShift(l, eta);
    const Real rhoStart = rho * farMultiplier;

    // WKB-corrected starting values at rhoStart, for a solution of the form
    // u(r) ~ sin(theta(r))/sqrt(k(r)) [regular-like, F] or
    // u(r) ~ cos(theta(r))/sqrt(k(r)) [irregular-like, G], theta(r) = the
    // standard Coulomb phase r - eta*ln(2r) - l*pi/2 + sigma_l(eta).
    // Amplitude AND its derivative both carry a 1/sqrt(k(r)) WKB
    // correction; using the leading-order (amplitude=1) form here was
    // empirically found, during implementation, to converge far too
    // slowly in farMultiplier to be practical -- see
    // docs/planning/coulomb-tail-continuum-matching.md.
    auto wkbStart = [&](Real r) {
        const Real theta = r - eta * std::log(2.0 * r) - l * kPi / 2.0 + sigmaL;
        const Real k = localK(r);
        const Real kp = localKPrime(r);
        const Real amp = 1.0 / std::sqrt(k);
        const Real sinT = std::sin(theta), cosT = std::cos(theta);
        const Real F0 = amp * sinT;
        const Real G0 = amp * cosT;
        const Real Fp0 = std::sqrt(k) * cosT - (kp / (2.0 * std::pow(k, 1.5))) * sinT;
        const Real Gp0 = -std::sqrt(k) * sinT - (kp / (2.0 * std::pow(k, 1.5))) * cosT;
        return std::array<Real, 4>{F0, Fp0, G0, Gp0};
    };

    // Exact Coulomb ODE, integrated inward (rhoStart > rho) via Numerov's
    // method: u'' = -coeff(r)*u, coeff(r) = 1 - 2*eta/r - l(l+1)/r^2 -- the
    // standard technique for a second-order ODE with no first-derivative
    // term (see the module-level comment on this function's declaration in
    // tise.hpp for why this replaced an initial 4th-order Runge-Kutta
    // attempt). Numerov's update needs two seed values (u at the two
    // farthest grid points, not value+derivative at one point), so the
    // second seed is obtained from a 3rd-order Taylor step off wkbStart's
    // value+derivative, using the ODE itself for the needed u''/u'''.
    auto coeffAt = [&](Real r) { return 1.0 - 2.0 * eta / r - l * (l + 1) / (r * r); };
    auto numerov = [&](Real u0, Real up0) {
        const Real diff = rhoStart - rho;
        const int nsteps = std::max(1, static_cast<int>(std::ceil(diff / stepSize)));
        const Real h = diff / nsteps; // exact fit: rhoStart - nsteps*h == rho
        const Real c0 = coeffAt(rhoStart);
        const Real u2 = -c0 * u0; // u'' = -coeff*u
        constexpr Real kEps = 1e-5;
        const Real cPrime = (coeffAt(rhoStart + kEps) - coeffAt(rhoStart - kEps)) / (2.0 * kEps);
        const Real u3 = -(cPrime * u0 + c0 * up0); // u''' = -(coeff'*u + coeff*u')
        Real rCurr = rhoStart - h;
        Real uPrev = u0; // at rhoStart (one point farther out than rCurr)
        Real uCurr = u0 - h * up0 + 0.5 * h * h * u2 - (h * h * h / 6.0) * u3; // Taylor step to rCurr
        Real fPrev = c0;
        Real fCurr = coeffAt(rCurr);
        for (int i = 0; i < nsteps - 1; ++i)
        {
            const Real rNext = rCurr - h;
            const Real fNext = coeffAt(rNext);
            const Real uNext = (2.0 * (1.0 - 5.0 * h * h * fCurr / 12.0) * uCurr -
                                 (1.0 + h * h * fPrev / 12.0) * uPrev) /
                                (1.0 + h * h * fNext / 12.0);
            uPrev = uCurr;
            uCurr = uNext;
            fPrev = fCurr;
            fCurr = fNext;
            rCurr = rNext;
        }
        // uCurr is now at rho (the target); uPrev is one step farther out, at
        // rho+h. Extract the derivative at rho from the Taylor expansion
        // uPrev = uCurr + h*u'(rho) + h^2/2*u''(rho), with u''(rho)=-fCurr*uCurr:
        //   u'(rho) = (uPrev - uCurr)/h + (h/2)*fCurr*uCurr
        const Real uPrimeCurr = (uPrev - uCurr) / h + 0.5 * h * fCurr * uCurr;
        return std::array<Real, 2>{uCurr, uPrimeCurr};
    };

    const auto start = wkbStart(rhoStart);
    const auto Fresult = numerov(start[0], start[1]);
    const auto Gresult = numerov(start[2], start[3]);

    return CoulombWaveResult{Fresult[0], Fresult[1], Gresult[0], Gresult[1]};
}

AsymptoticResult matchAsymptotic(
    const bspline::BSpline &bs,
    std::vector<std::vector<Real>> states,
    const EigenResult &eigen,
    std::vector<Real> grid, Real R,
    int order, const std::vector<Real> &Hmat, const std::vector<Real> &Smat,
    std::optional<std::vector<int>> dropSet,
    Real fineDE,
    std::optional<std::pair<int, Real>> coulombLC
)
{
    AsymptoticResult result;
    result.A_E = std::vector<Real>(grid.size(), 0.0);
    result.delta = std::vector<Real>(grid.size(), 0.0);
    result.dDeltaDE = std::vector<Real>(grid.size(), 0.0);

    int nEn = eigen.dim;
    int nBSplines = bs.getNBSplines();

    // Generalized via colOf/physicalOf (A4b, REQ-F-050) instead of the
    // classic hardcoded fc[1+j]/fc[nEn+1] shift-by-1 placement: that offset
    // only works when exactly B_1 is dropped (every kept physical index i
    // sits at coefficient slot i-1). With interior B-splines also droppable,
    // physicalOf[j] gives the true physical index kept eigenvector column j
    // corresponds to, so continuumStateToBSplineCoeffs can place each
    // coefficient at the right physical B-spline regardless of which
    // interior indices were dropped.
    const ColumnMap map = columnIndexMap(nBSplines, dropSet.value_or(std::vector<int>{1}));
    const std::vector<int> &physicalOf = map.physicalOf;

    // Needed only for the fine-step dDeltaDE evaluations below (rawDeltaAt)
    // -- cheap (a boundary-coupling solve, not a re-diagonalization), and
    // independent of the passed-in `states` (built only at the production
    // grid, which the fine steps fall outside of).
    auto [coeffs1, coeffs2] = precomputeBoundaryCoupling(order, nEn, Hmat, Smat, eigen, nBSplines, dropSet);

    // Value+derivative match at R against A_E*sin(kR+delta) (flat) or
    // A_E[cos(delta)*F_l(eta,kR)+sin(delta)*G_l(eta,kR)] (Coulomb, when
    // coulombLC is set) -- see the module-level comment above
    // matchAsymptotic's declaration in tise.hpp. Both branches solve for
    // the SAME quantities (A_E, delta) via the value/derivative match;
    // the Coulomb branch's Wronskian-based solve reduces to the flat
    // branch's atan/sqrt formulas exactly at eta=0 (F_0(0,rho)=sin(rho),
    // G_0(0,rho)=cos(rho)), confirmed during implementation.
    //
    // eta = C/k is recomputed fresh from k EVERY call (not hoisted/fixed
    // once outside this lambda) -- eta genuinely depends on the energy via
    // k=sqrt(2E), which varies across the energy grid this lambda is
    // called once per point of.
    auto amplitudeAndDelta = [&](Real psi_R, Real psiPrime_R, Real k) -> std::pair<Real, Real> {
        if (!coulombLC)
        {
            const Real A_E = std::sqrt((2.0 / M_PI) / (k * psi_R * psi_R + psiPrime_R * psiPrime_R / k));
            const Real delta = std::atan(k * psi_R / psiPrime_R) - k * R;
            return {A_E, delta};
        }
        const auto [l, C] = *coulombLC;
        const Real eta = C / k;
        const CoulombWaveResult cw = evaluateCoulombFunctions(l, eta, k * R);
        const Real psiPrimeOverK = psiPrime_R / k;
        const Real W = cw.F * cw.Gprime - cw.G * cw.Fprime; // ~1 by construction; computed, not assumed
        const Real alpha = (psi_R * cw.Gprime - psiPrimeOverK * cw.G) / W; // A_E*cos(delta)
        const Real beta = (psiPrimeOverK * cw.F - psi_R * cw.Fprime) / W; // A_E*sin(delta)
        const Real A_E = std::sqrt((2.0 / (M_PI * k)) / (alpha * alpha + beta * beta));
        const Real delta = std::atan2(beta, alpha);
        return {A_E, delta};
    };

    // for each E, first find \bar psi_E(R) and \bar psi'_E(R), then calculate A_E, delta
    for (int E_idx = 0; E_idx < grid.size(); ++E_idx)
    {
        // states[E_idx][i] (i=0..nEn-1) give the coefficients for bound eigenstate phi_i which builds |\bar psi_E>,
        // (and states[E_idx][nEn] is always 1.0, the coefficient for B_N)
        // To compute the scalar product needed to find \bar psi_E(R) and \bar psi'_E(R) with bs.eval, however,
        // we need ALL coefficients to correspond to the B-Splines, not phi_n!!!
        std::vector<Real> fc = continuumStateToBSplineCoeffs(states[E_idx], eigen, nBSplines, physicalOf);

        // NOW we can use bs.eval
        Real psi_R = bs.eval(R, fc.data(), fc.size(), 0);
        Real psiPrime_R = bs.eval(R, fc.data(), fc.size(), 1);

        Real k = sqrt(2 * grid[E_idx]);

        std::tie(result.A_E[E_idx], result.delta[E_idx]) = amplitudeAndDelta(psi_R, psiPrime_R, k);
    }

    // === H2/H3 fix (docs/tests/reports/8236239/finite_square_well.md,
    // hydrogen.md): the previous dDeltaDE finite-differenced sin(2*delta)
    // across the PRODUCTION energy grid and divided by 2*cos(2*delta). That
    // failed two independent ways: the production grid is far too coarse
    // for delta's actual variation (routinely O(1) rad between points), and
    // the sin/cos inversion is an inherent 0/0 conditioning problem at
    // delta=pi/4 (mod pi/2) that no step size fixes. Replaced by: (1)
    // storing a CONTINUOUS, unwrapped delta(E) instead of raw
    // atan(...)-kR, and (2) computing dDeltaDE as a central difference of
    // that same raw delta over a small internal step around each requested
    // energy, independent of the production grid's own spacing.
    auto rawDeltaAt = [&](Real E) -> Real {
        std::vector<Real> s(nEn + 1, 0.0);
        for (int i = 0; i < nEn; ++i)
            s[i] = (coeffs1[i] - E * coeffs2[i]) / (E - eigen.values[i]);
        s[nEn] = 1.0;
        std::vector<Real> fc = continuumStateToBSplineCoeffs(s, eigen, nBSplines, physicalOf);
        Real psi_R = bs.eval(R, fc.data(), fc.size(), 0);
        Real psiPrime_R = bs.eval(R, fc.data(), fc.size(), 1);
        Real k = std::sqrt(2 * E);
        return amplitudeAndDelta(psi_R, psiPrime_R, k).second;
    };

    for (std::size_t E_idx = 0; E_idx < grid.size(); ++E_idx)
    {
        const Real E = grid[E_idx];
        const Real h = std::min(fineDE, 0.5 * E); // keep E-h strictly > 0
        Real deltaMinus = rawDeltaAt(E - h);
        Real deltaPlus  = rawDeltaAt(E + h);
        // Unwrap this tiny window relative to itself (atan's branch jumps
        // by pi when psiPrime_R crosses zero) -- a resonance landing
        // exactly inside so small a window is already surfaced separately
        // by buildContinuumState's pole-proximity warning.
        Real diff = deltaPlus - deltaMinus;
        while (diff > M_PI / 2)  { deltaPlus -= M_PI; diff -= M_PI; }
        while (diff < -M_PI / 2) { deltaPlus += M_PI; diff += M_PI; }
        result.dDeltaDE[E_idx] = (deltaPlus - deltaMinus) / (2.0 * h);
    }

    // Unwrap the stored production-grid delta(E) itself into a continuous
    // trajectory (closes H3 -- docs/tests/reports/8236239/free_particle.md:
    // raw atan(...)-kR runs to -100+ and makes phase_shifts.png
    // "uninformative as drawn"). Normalize the first point into
    // (-pi/2, pi/2], then track pi-periodic branch jumps thereafter --
    // using a derivative-informed (not just closest-to-previous-raw-value)
    // prediction: delta genuinely changes by more than pi/2 between
    // production grid points for some potentials (e.g. finite_square_well
    // at its default n_energies:5 grid, dDeltaDE ~ -13 to -19 per unit E,
    // ~0.1 apart -- an actual change of ~1.5-2 rad, comparable to or larger
    // than the pi/2 "closest branch" heuristic can disambiguate on its
    // own). result.dDeltaDE (computed above, from the fine-step method,
    // independent of production grid spacing) gives an accurate LOCAL
    // slope at each already-unwrapped point; a first-order Euler
    // extrapolation from it predicts the next point far more reliably than
    // assuming the true change is small.
    if (!result.delta.empty())
    {
        while (result.delta[0] > M_PI / 2)   result.delta[0] -= M_PI;
        while (result.delta[0] <= -M_PI / 2) result.delta[0] += M_PI;
        for (std::size_t E_idx = 1; E_idx < result.delta.size(); ++E_idx)
        {
            const Real predicted = result.delta[E_idx - 1] +
                result.dDeltaDE[E_idx - 1] * (grid[E_idx] - grid[E_idx - 1]);
            Real diff = result.delta[E_idx] - predicted;
            while (diff > M_PI / 2)  { result.delta[E_idx] -= M_PI; diff -= M_PI; }
            while (diff < -M_PI / 2) { result.delta[E_idx] += M_PI; diff += M_PI; }
        }
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
                     Real rMax,
                     const EigenResult &eigen,
                     std::optional<std::vector<int>> dropSet)
{
    // phase_shifts.dat: 3-col epsilon_i, delta(epsilon_i), dDelta/dE
    out << std::scientific << std::setprecision(16);
    for (std::size_t i = 0; i < grid.size(); ++i)
    {
        out << " " << std::setw(24) << grid[i]
            << " " << std::setw(24) << result.delta[i]
            << " " << std::setw(24) << result.dDeltaDE[i] << "\n";
    }

    const int nBSplines = bs.getNBSplines();
    const ColumnMap map = columnIndexMap(nBSplines, dropSet.value_or(std::vector<int>{1}));
    const std::vector<int> &physicalOf = map.physicalOf;

    // continuum_state_NNN.dat: 2-col x, psi_{epsilon_i}(x), one block per energy
    for (std::size_t i = 0; i < grid.size(); ++i)
    {
        // states[i] holds coefficients of the confined eigenstates {phi_n}
        // (plus B_N) -- transform into true B-spline coefficients before
        // calling bs.eval, exactly like matchAsymptotic does (same shared
        // helper, not a second copy of the transform).
        std::vector<Real> fc = continuumStateToBSplineCoeffs(states[i], eigen, nBSplines, physicalOf);

        std::ostream &stOut = *stateOut[i];
        stOut << std::scientific << std::setprecision(16);
        for (int ix = 1; ix <= npts; ++ix)
        {
            Real x = rMin + (rMax - rMin) *
                     static_cast<Real>(ix - 1) / static_cast<Real>(npts - 1);
            Real psi = result.A_E[i] * bs.eval(x, fc.data(), fc.size(), 0);
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
// instead of growing toward infinity. This is valid reuse because
// classifySequenceConvergence's divergence test only looks at how finite
// differences shrink/grow across the array it's given -- it has no built-in
// assumption about whether the sample points are marching outward toward
// infinity (classifyAsymptote's use) or inward toward a finite point (this
// function's use); either way, a divergent potential produces
// non-shrinking successive differences. `sign` is -1 to sample the piece
// left of x0, +1 for the piece right of x0; `pieceWidth` is that piece's
// own width, used to clamp the starting offset so the probe can never step
// outside a narrow piece (e.g. a thin barrier slab) -- classifyAsymptote's
// own scale formula alone isn't safe here since it assumes room to grow
// outward, not a bounded piece to stay inside of.
//
// numSamples/ratio/backupScale share classifyAsymptote's defaults (same
// probing idea, just inverted direction); pieceWidthClampFraction is this
// function's own knob on how much of the piece's width the first probe step
// may use.
bool isSingularApproaching(const std::map<std::string, std::string> &potential,
                            Real x0, Real sign, Real pieceWidth,
                            int numSamples = kDefaultAsymptoteNumSamples,
                            Real ratio = kDefaultAsymptoteRatio,
                            Real backupScale = kDefaultAsymptoteBackupScale,
                            Real pieceWidthClampFraction = 0.25)
{
    // Same requirement as classifyAsymptote's own guard: numSamples becomes
    // V's size below, and classifySequenceConvergence's tail-window fit
    // needs at least 3 samples to produce one estimate -- validated here,
    // before any sampling work, for the same fail-fast/clear-error reason.
    if (numSamples < 3)
        throw std::runtime_error("isSingularApproaching: numSamples must be >= 3");

    const Real scale = std::max(std::abs(x0), backupScale);
    // Never let the first probe step exceed pieceWidthClampFraction of the
    // piece's own width, so a narrow piece can't have its probe overshoot
    // into (or past) the neighboring piece.
    const Real h0 = std::isfinite(pieceWidth) ? std::min(scale, pieceWidthClampFraction * pieceWidth) : scale;

    std::vector<Real> V(numSamples);
    for (int k = 0; k < numSamples; ++k)
        V[k] = evaluateFunction(potential, x0 + sign * h0 / std::pow(ratio, k));

    return classifySequenceConvergence(V, ratio).isDivergent;
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

std::vector<DetectedJoin> detectPotentialStructure(const std::map<std::string, std::string> &potential,
                                                     Real valueJumpTol,
                                                     Real slopeJumpTol,
                                                     Real fdStep,
                                                     Real stepClampFraction)
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

    // valueJumpTol: O(1)-O(10) potential magnitudes in this codebase; a
    // genuine Step is an O(1) jump, so the default (1e-6) leaves ~6 orders
    // of margin against finite-difference noise.
    // slopeJumpTol: the one-sided slope estimate's own O(h^2) truncation
    // error is ~1e-8 at the default fdStep=1e-4, ~4 orders below the
    // default slopeJumpTol (1e-4).
    // stepClampFraction: never let the finite-difference step exceed this
    // fraction of a piece's own width, so probing near a join can't step
    // outside a narrow piece into its neighbor.
    auto clampedStep = [&](Real width) {
        return std::isfinite(width) ? std::min(fdStep, stepClampFraction * width) : fdStep;
    };

    auto classifyJoin = [&](Real x0, Real leftWidth, Real rightWidth) {
        if (isSingularApproaching(potential, x0, -1.0, leftWidth) ||
            isSingularApproaching(potential, x0, +1.0, rightWidth))
            return JoinType::Singular;

        const Real hL = clampedStep(leftWidth);
        const Real hR = clampedStep(rightWidth);
        const Real vL = oneSidedValue(potential, x0, -1.0, hL);
        const Real vR = oneSidedValue(potential, x0, +1.0, hR);
        if (std::abs(vL - vR) > valueJumpTol)
            return JoinType::Step;

        const Real sL = oneSidedSlope(potential, x0, -1.0, hL);
        const Real sR = oneSidedSlope(potential, x0, +1.0, hR);
        if (std::abs(sL - sR) > slopeJumpTol)
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

// === A4b: generalized B-spline drop-set (REQ-F-050 table row 4) ===
// Physical B-spline indices (1-based) whose support touches x, given the
// same (nNodes, order, grid) that will be passed to BSpline::init. Support
// of B-spline Bs spans extended-knot indices [Bs-order, Bs] (verified
// against BSpline::init's own construction, BSpline.cpp:204-230: the
// extended grid clamps to grid.front()/grid.back() outside [0,nNodes-1],
// and each B-spline's defining knot vector starts at extended index
// Bs-order). Closed-interval ("touches", not strict interior support): a
// B-spline whose support ends exactly at x is included -- appropriate for
// identifying removal candidates at exactly a singular x, where being
// slightly inclusive is the conservative choice. Note: touching a domain
// *boundary* point pulls in an entire cluster of `order` B-splines (all of
// them clamp to the same boundary knot), not just one. Worked example
// (order=4, x = grid.front()): the extended knot vector repeats the
// boundary knot `order` times, so B-splines 1..order all have support
// beginning exactly at x and every one of them "touches" it by this
// closed-interval test, not just B_1.
std::vector<int> bSplinesTouchingX(int nNodes, int order, const std::vector<Real> &grid, Real x)
{
    auto knotAt = [&](int i) -> Real {
        if (i < 0) return grid.front();
        if (i >= nNodes) return grid.back();
        return grid[i];
    };

    const int nBSplines = nNodes + order - 2;
    std::vector<int> touching;
    for (int bs = 1; bs <= nBSplines; ++bs)
    {
        const Real lo = knotAt(bs - order);
        const Real hi = knotAt(bs);
        if (x >= lo && x <= hi)
            touching.push_back(bs);
    }
    return touching;
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
                                          int nBSplines,
                                          std::optional<std::vector<int>> dropSet)
{
    const ColumnMap map = resolveDropSet(nBSplines, nEn, dropSet);
    const std::vector<int> &colOf = map.colOf;

    std::vector<Real> coeffs(nBSplines, 0.0);
    const Real *col = &evec[(iEn - 1) * nEn];
    for (int physicalIdx = 1; physicalIdx <= nBSplines; ++physicalIdx)
        if (colOf[physicalIdx] != 0)
            coeffs[physicalIdx - 1] = col[colOf[physicalIdx] - 1];
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

void writeEigenvalues(std::ostream &out, const EigenResult &er, int nStates)
{
    out << "# eigenvalues.dat: index, E_n\n";
    out << std::scientific << std::setprecision(16);
    for (int i = 0; i < nStates; ++i)
        out << " " << std::setw(6) << i
            << " " << std::setw(24) << er.values[i] << "\n";
}

void writeEigenvectors(std::ostream &out, const EigenResult &er, int nBSplines, int nStates,
                        std::optional<std::vector<int>> dropSet)
{
    std::vector<std::vector<Real>> cols(nStates);
    for (int j = 0; j < nStates; ++j)
        cols[j] = eigenstateCoefficients(er.vectors, j + 1, er.dim, nBSplines, dropSet);

    out << "# eigenvectors.dat: columns are c_n coefficient vectors\n";
    out << std::scientific << std::setprecision(16);
    for (int row = 0; row < nBSplines; ++row)
    {
        for (int j = 0; j < nStates; ++j)
            out << " " << std::setw(24) << cols[j][row];
        out << "\n";
    }
}

void writeBandedMatrix(std::ostream &out, const std::vector<Real> &mat,
                        int order, int nEn, const std::string &description)
{
    out << "# " << description << "\n";
    out << std::scientific << std::setprecision(16);
    for (int row = 0; row < order; ++row)
    {
        for (int col = 0; col < nEn; ++col)
            out << " " << std::setw(24) << mat[row + col * order];
        out << "\n";
    }
}

StrategicGridResult buildStrategicGridAndDropSet(int nNodes, int order, Real rMin, Real rMax,
                                                   const std::map<std::string, std::string> &potential,
                                                   Real edgeTolerance)
{
    // Strategic grid construction is automatic, not something a caller
    // opts into separately: detectPotentialStructure always runs first, and
    // its result alone decides whether the grid ends up strategic. A
    // potential with no detectable Step/StitchedKink structure produces an
    // unchanged uniform grid, byte-identical to a plain buildUniformRadialGrid
    // call.
    auto joins = detectPotentialStructure(potential);
    auto knots = strategicKnotsFromJoins(joins, order);

    // docs/tests/reports/8236239/interior_singularity.md: a genuine INTERIOR
    // Singular join (not already regularized by a domain edge) needs knot
    // multiplicity order-1 -- one short of a full clamp -- so that exactly
    // ONE B-spline is non-zero there, mirroring how a domain edge's clamped
    // knots leave exactly one boundary B-spline (B_1/B_N) non-zero. Without
    // this, the join sits on an ordinary simple knot and the fillDropSet
    // loop below (bSplinesTouchingX) would drop the entire ~order-sized
    // touching cluster instead of the single point value, forcing every
    // eigenfunction to zero over a multi-bohr "dead zone" around the join --
    // confirmed to produce eigenvalues 41% too high on the field-free side
    // of a box/repulsive-Coulomb split. Edge-Singular joins are excluded
    // here (extraMultiplicity 0, same as strategicKnotsFromJoins already
    // gives every Singular join) -- they're untouched, exactly as before.
    for (const auto &j : joins)
    {
        if (j.type != JoinType::Singular)
            continue;
        const bool atLeftEdge  = std::abs(j.x - rMin) < 1e-9;
        const bool atRightEdge = std::abs(j.x - rMax) < 1e-9;
        if (atLeftEdge || atRightEdge)
            continue;
        if (order - 2 > 0)
            knots.push_back({j.x, order - 2});
    }

    auto grid  = buildStrategicRadialGrid(nNodes, rMin, rMax, knots);
    const int nNodesActual = static_cast<int>(grid.size());

    bspline::BSpline bs;
    int initInfo = bs.init(nNodesActual, order, grid);
    if (initInfo != 0)
        throw std::runtime_error("BSpline::init failed with code " +
                                 std::to_string(initInfo));

    int nBSplines = bs.getNBSplines();

    // Gap 2: singular-join B-spline removal -- INTERIOR joins only. A
    // Singular join located AT a domain edge (e.g. hydrogen's Coulomb
    // singularity coinciding with the left wall at x=rMin) is already
    // regularized by the classic single-B-spline wall exclusion (that's
    // exactly what "drop B_1" already enforces: u(rMin)=0). Removing a
    // FULL bSplinesTouchingX cluster there instead -- verified empirically
    // during implementation -- guts the basis precisely where a
    // near-origin-peaked wavefunction (e.g. hydrogen's ground state) has
    // most of its amplitude, producing a qualitatively wrong result (a 5x
    // ground-state error was observed for a boundary-cluster removal at
    // order=8), not just reduced accuracy. A genuine INTERIOR singularity
    // (no piece of the domain boundary already regularizes it) is the case
    // this mechanism is actually for.
    //
    // For an interior join, dropping the FULL touching cluster (as this
    // used to do) is likewise too aggressive: on the simple-knot grid that
    // used to be built here, ~order B-splines touch the point, and removing
    // all of them forces every eigenfunction to zero over a multi-bohr dead
    // zone around it (docs/tests/reports/8236239/interior_singularity.md,
    // confirmed 41% eigenvalue error). The knot-multiplicity bump added
    // above leaves exactly ONE B-spline non-zero at the join -- find it
    // numerically (bs.eval at j.x for each candidate the multiplicity-aware
    // grid now makes bSplinesTouchingX return) and drop only that one,
    // exactly mirroring the single-B-spline edge treatment.
    //
    // B_N (nBSplines) is deliberately never added to fillDropSet even if a
    // right-edge cluster would include it -- continuum construction needs
    // its raw H/S column filled, not dropped; it's excluded from the
    // DIAGONALIZATION instead, by truncation below (generalizes the
    // pre-existing {1}/nEn+1/truncate-to-nEn pattern to an arbitrary-size
    // interior cluster).
    std::vector<int> fillDropSet = {1};
    bool rightEdgeSingular = false;
    for (const auto &j : joins)
    {
        if (j.type != JoinType::Singular)
            continue;

        const bool atLeftEdge  = std::abs(j.x - rMin) < edgeTolerance;
        const bool atRightEdge = std::abs(j.x - rMax) < edgeTolerance;

        // matchAsymptotic's flat-asymptote assumption at R~rMax is affected
        // by a right-edge singularity regardless of whether any extra
        // B-spline removal happens there -- always flag it for the caller.
        if (atRightEdge)
            rightEdgeSingular = true;

        if (atLeftEdge || atRightEdge)
            continue; // already regularized by the classic wall exclusion

        auto candidates = bSplinesTouchingX(nNodesActual, order, grid, j.x);
        constexpr Real kNonzeroTol = 1e-6; // B-spline peak values are O(1);
                                            // roundoff for a mathematically-
                                            // zero candidate is ~1e-12.
        for (int idx : candidates)
            if (std::abs(bs.eval(j.x, idx, 0)) > kNonzeroTol)
                fillDropSet.push_back(idx);
    }
    std::sort(fillDropSet.begin(), fillDropSet.end());
    fillDropSet.erase(std::unique(fillDropSet.begin(), fillDropSet.end()), fillDropSet.end());
    fillDropSet.erase(std::remove(fillDropSet.begin(), fillDropSet.end(), nBSplines), fillDropSet.end());

    const int nEnFilled = nBSplines - static_cast<int>(fillDropSet.size()); // == nEnBound + 1
    const int nEnBound  = nEnFilled - 1;

    return StrategicGridResult{grid, bs, nBSplines, fillDropSet, nEnBound, rightEdgeSingular};
}

SolveTISEResult solveTISE(int nNodes, int order, Real rMin, Real rMax, int L, std::map<std::string, std::string> potential,
                           Real E_threshold, Real E_max, int N_E, int continuumOutputPoints)
{
    auto sgr = buildStrategicGridAndDropSet(nNodes, order, rMin, rMax, potential);

    if (sgr.rightEdgeSingular)
        std::cerr << "Warning: potential is singular at the right domain edge x=" << rMax
                  << "; continuum phase-shift matching (matchAsymptotic) assumes a regular "
                  << "boundary there, so continuum results should be treated with suspicion.\n";

    const int nEnFilled = sgr.nEnBound + 1;

    auto [H, S] = fillBandedMatrices(sgr.bs, nEnFilled, order, L, potential, sgr.fillDropSet);
    EigenResult er = solveGeneralizedEigenproblem(H, S, sgr.nEnBound, order);

    auto energyGrid = buildEnergyGrid(E_threshold, E_max, N_E);
    std::vector<std::vector<tise::Real>> states =
        buildContinuumState(order, sgr.nEnBound, H, S, er, energyGrid, sgr.nBSplines, sgr.fillDropSet);
    AsymptoticResult ar = matchAsymptotic(sgr.bs, states, er, energyGrid, rMax, order, H, S, sgr.fillDropSet);

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
    writeContinuumInfo(phaseShiftsOut, sgr.bs, ar, energyGrid, states, continuumStateOut,
                        continuumOutputPoints, rMin, rMax, er, sgr.fillDropSet);

    std::vector<int> dropSet = sgr.fillDropSet;
    dropSet.push_back(sgr.nBSplines);
    std::sort(dropSet.begin(), dropSet.end());

    return SolveTISEResult{er, sgr.bs, sgr.grid, sgr.nBSplines, dropSet};
}

// === A5: E_acc continuum-accuracy warning (REQ-F-040, warning half) ===
Real minInterNodeGap(const std::vector<Real> &grid, Real tol)
{
    if (grid.size() < 2)
        throw std::runtime_error("minInterNodeGap: grid must have at least 2 points");
    Real minGap = std::numeric_limits<Real>::infinity();
    Real lastDistinct = grid[0];
    for (std::size_t i = 1; i < grid.size(); ++i)
    {
        const Real gap = grid[i] - lastDistinct;
        if (gap > tol)
        {
            minGap = std::min(minGap, gap);
            lastDistinct = grid[i];
        }
    }
    if (!std::isfinite(minGap))
        throw std::runtime_error("minInterNodeGap: fewer than 2 distinct points in grid");
    return minGap;
}

Real computeEAcc(Real nodeSpacing, Real mass)
{
    // E = pi^2 / (2 * mass * nodeSpacing^2), straight from the derivation in
    // this function's own tise.hpp doc comment (k >~ pi/dx_node, E = k^2/2m).
    // The `2.0` here is the standard kinetic-energy mass factor from
    // E = p^2/(2m) = (hbar*k)^2/(2m) in atomic units (hbar=1) -- a fixed
    // term of the physics formula itself, not a tunable numerical knob, so
    // it stays a literal rather than becoming a named constant/parameter.
    return kPi * kPi / (2.0 * mass * nodeSpacing * nodeSpacing);
}

bool warnIfContinuumExceedsEAcc(Real eMax, Real eAcc, std::ostream &warnOut)
{
    const bool exceeds = eMax > eAcc;
    if (exceeds)
        warnOut << "Warning: requested continuum E_max=" << eMax
                << " exceeds the basis accuracy ceiling E_acc=" << eAcc
                << " (set by the B-spline node spacing); results at energies "
                << "above E_acc are unreliable.\n";
    return exceeds;
}

} // namespace tise
