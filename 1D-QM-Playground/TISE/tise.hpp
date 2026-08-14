#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <map>
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

// True if x falls within `interval`, e.g. "[0, 20)", "(-inf, inf)". Bounds may
// be finite numbers or (+/-)inf/infinity (case-insensitive); brackets select
// inclusive ([, ]) vs. exclusive ((, )) endpoints. Throws std::runtime_error
// if `interval` doesn't match the expected syntax.
bool inInterval(double x, const std::string &interval);

// Given a piecewise potential (domain string -> muparser expression string in
// `x`), find the piece whose domain contains x and evaluate it there. Throws
// std::runtime_error if no piece's domain covers x.
double evaluateFunction(std::map<std::string, std::string> function, double x);

// Result of fitting a sampled sequence V[0..N-1] (assumed to be sampled at
// points growing/shrinking geometrically by `ratio` per step) for convergence
// behavior: does it diverge, is it already flat, or does it follow a
// power-law V ~ fittedLimit + C/step^powerLawExponent?
struct ConvergenceFit
{
    bool isDivergent;
    bool isFlat;
    Real fittedLimit;
    Real powerLawExponent;
};

// Shared numeric core for asymptote classification (REQ-F-030) and, later,
// singular-point detection for strategic node placement (REQ-F-050). See
// docs/planning/engineer-a-plan-A1.md for the algorithm.
ConvergenceFit classifySequenceConvergence(const std::vector<Real> &V, Real ratio);

// Which side of the spatial domain a boundary/asymptote calculation refers to.
enum class DomainSide
{
    Left,
    Right
};

// Smooth sin^2 raised-cosine taper used for Case-3 (irregular-asymptote)
// boundary handling, per docs/planning/boundary-condition-case-3-smoothing.md:
// 1 in the trusted interior, 0 at and beyond the boundary R, C^1-continuous
// in between over a transition width `delta`. For DomainSide::Right, the
// interior is x <= R-delta; for DomainSide::Left, the interior is x >= R+delta.
Real case3WindowFunction(Real x, Real R, Real delta, DomainSide side);

// V~(x) = W(x) * V(x): the smoothed, boundary-tapered potential used in place
// of a hard flat-truncation for Case 3. `potential` is evaluated via
// evaluateFunction; the result is then multiplied by case3WindowFunction.
Real evaluateWindowedPotential(const std::map<std::string, std::string> &potential,
                                Real x, Real R, Real delta, DomainSide side);

// The spatial domain of the TISE problem (a single finite box; Dirichlet
// walls at both xMin and xMax regardless of which sides are "unbounded" in
// the physical problem being approximated).
struct SpatialDomain
{
    Real xMin;
    Real xMax;
};

// Figure 7 (SDD Sec. 5.2.3) Case 1/2/3 dispatch for an unbounded domain side.
enum class AsymptoteCase
{
    HardWall,          // Case 1: no finite asymptote (diverges/grows)
    AnalyticAsymptote, // Case 2: flat or Coulomb ~1/r
    Irregular          // Case 3: unknown/irregular tail
};

// Case-2 sub-branch. Only Flat's continuum-matching formula is implemented
// (Engineer B's B3); Coulomb is classified but its matching formula is not
// yet derived by any source document (flagged as a follow-up).
enum class AsymptoteSubType
{
    NotApplicable,
    Flat,
    Coulomb
};

struct AsymptoteClassification
{
    AsymptoteCase    asymptoteCase;
    AsymptoteSubType subType;
    Real             fittedAsymptoticValue;      // V_inf; NaN if HardWall
    Real             powerLawExponent;           // fitted p; NaN if HardWall or Flat
    Real             recommendedTransitionWidth; // Delta for the Case-3 window; NaN unless Irregular
    bool             warningEmitted;
};

// Classify the potential's asymptote on the given unbounded domain side
// (REQ-F-030). Assumes the caller already knows this side is unbounded;
// bounded sides always get a plain Dirichlet wall and never call this.
// Case-3 warnings are written to `warnOut` (physics warning, non-fatal,
// per SDD Sec. 8's warning taxonomy).
AsymptoteClassification classifyAsymptote(const std::map<std::string, std::string> &potential,
                                           const SpatialDomain &domain,
                                           DomainSide side,
                                           std::ostream &warnOut = std::cerr);

// Fill symmetric banded Hamiltonian H and overlap S matrices (LAPACK 'U' storage).
// Returns {Hmat, Smat}, each of size order * nEn.
// `potential` maps domain strings (e.g. "[0, inf)") to muparser expressions in
// `x`; the piece whose domain contains a given x is evaluated to give V(x).
std::pair<std::vector<Real>, std::vector<Real>>
fillBandedMatrices(const bspline::BSpline &bs, int nEn, int order, int L, std::map<std::string, std::string> potential);

// Solve H c = E S c via LAPACK DSBGV.
// H and S are consumed (overwritten); pass by value intentionally.
EigenResult solveGeneralizedEigenproblem(std::vector<Real> H,
                                         std::vector<Real> S,
                                         int nEn,
                                         int order);

// === A2: bound-state classification (REQ-F-020) ===
// Classification of a solved eigenproblem's eigenvalues against an
// ionization threshold supplied by the caller (the threshold itself is
// never computed here -- see classifyAsymptote for that, or a caller may
// know it by convention, e.g. 0.0 for a Coulomb-type problem). Does not
// mutate EigenResult.
struct BoundStateClassification
{
    std::vector<bool> isBound; // isBound[i] true iff values[i] < threshold
    int nBound;                // count of true entries in isBound
};

// Strict less-than: a state exactly at threshold is the marginal case and
// is classified as an above-threshold pseudostate, not bound.
BoundStateClassification classifyBoundStates(const EigenResult &result, Real threshold);

// === A3: well-containment diagnostic (SDD Sec. 5.2.3, Sec. 6.4) ===
// Result of checkWellContainment: the raw signed derivative plus the
// pass/fail flag derived from it.
struct ContainmentCheck
{
    Real psiPrimeAtBoundary; // psi'(xBoundary), as returned by BSpline::eval
    bool notWellContained;   // true iff |psiPrimeAtBoundary| > tol
};

// Well-containment diagnostic (SDD Sec. 5.2.3 Figure 6 "BOUND" node; Sec. 6.4
// data-validation rule): a bound state confined within the box should decay
// to numerically-zero slope at the outer wall; a nonzero psi'(xBoundary)
// means the state is "colliding" with the wall and its energy/wavefunction
// may be inaccurate. `coeffs` is whatever eigenstateCoefficients produces
// (already zero-padded at both ends) for a single eigenstate; `xBoundary` is
// typically rMax (or rMin for a left-side check). Strict '>': a derivative
// magnitude exactly equal to tol is the marginal case and is not flagged,
// mirroring classifyBoundStates' strict-inequality boundary convention.
ContainmentCheck checkWellContainment(const bspline::BSpline &bs,
                                       const std::vector<Real> &coeffs,
                                       Real xBoundary,
                                       Real tol = 1e-3);

// === A4a: strategic node placement (REQ-F-050), required core only ===
// A4b (generalizing the hardcoded "drop B_1/B_N" convention so singular-
// potential B-spline removal is more than detection-only) is deferred --
// see docs/planning/engineer-a-plan.md Section A4. Singular joins are
// detected here but not remediated by this task.

// One of the four required-treatment categories from the task's source
// table (docs/planning/tise-task-breakdown.md, "A4. Strategic node
// placement"): Step (V itself jumps), StitchedKink (V continuous, V' jumps),
// Singular (V diverges, e.g. 1/r), or Continuous (no special treatment).
enum class JoinType
{
    Continuous,
    Step,
    StitchedKink,
    Singular
};

// A detected join (interior piece-to-piece boundary) or flagged domain edge
// in the potential's piecewise structure, at position x.
struct DetectedJoin
{
    Real x;
    JoinType type;
};

// A single knot location that should receive `extraMultiplicity` additional
// degenerate copies beyond the ordinary grid point already there (or
// spliced in if not already present -- see buildStrategicRadialGrid).
struct StrategicKnot
{
    Real x;
    int extraMultiplicity;
};

// Detect step/stitched-kink/singular structure in a piecewise potential's
// domain joins and finite outer domain edges. Interior joins (two pieces
// meet) are always emitted, even Continuous ones -- so a Continuous join's
// presence in the output means "detected, not actionable", not "silently
// dropped". Finite domain edges (a piece's own outer boundary with no
// neighbor, e.g. x=0 in {"(0,100]": "-1/x+1/x^2"}) are only emitted if
// singular there -- an ordinary box wall isn't actionable for this task.
// See docs/planning/engineer-a-plan-A4.md for the detection algorithm and
// threshold derivations.
std::vector<DetectedJoin> detectPotentialStructure(const std::map<std::string, std::string> &potential);

// Convert detected joins into strategic knots, per the multiplicity formula
// derived in docs/planning/engineer-a-plan-A4.md: Step -> order-3,
// StitchedKink -> order-4 (both clamped at 0 -- a case this codebase's
// order values, 4 and 8, only actually reach for StitchedKink at order=4).
// Singular/Continuous produce no knot (Singular's remediation -- B-spline
// removal -- is A4b, out of scope here).
std::vector<StrategicKnot> strategicKnotsFromJoins(const std::vector<DetectedJoin> &joins, int order);

// Build a radial grid combining a uniform base (buildUniformRadialGrid) with
// degenerate knots at each StrategicKnot's location: extraMultiplicity
// repeated copies, splicing in a brand-new grid point first if `x` doesn't
// already coincide (within 1e-12) with a uniform grid point. Returns the
// FULL combined grid.
//
// CONTRACT: callers MUST pass the returned vector's own .size() as
// BSpline::init's `numberOfNodes` argument, not the original nNodes --
// BSpline::init only reads the first numberOfNodes entries of the grid it's
// given and silently drops the rest (verified against BSpline.cpp:190-216).
std::vector<Real> buildStrategicRadialGrid(int nNodes, Real rMin, Real rMax,
                                            const std::vector<StrategicKnot> &knots);

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
// `potential` is the piecewise potential passed through to fillBandedMatrices.
EigenResult solveTISE(int nNodes, int order, Real rMin, Real rMax, int L, std::map<std::string, std::string> potential);

} // namespace tise
