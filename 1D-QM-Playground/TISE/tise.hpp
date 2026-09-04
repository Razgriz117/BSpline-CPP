#pragma once

#include <iostream>
#include <optional>
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

// Holds the output of matchAsymptotic, A_E, delta, dDeltaDE for each energy E in the input energy grid.
struct AsymptoticResult
{
    std::vector<Real> A_E;
    std::vector<Real> delta;
    std::vector<Real> dDeltaDE;
};

// Build a uniform radial grid of nNodes points on [rMin, rMax].
std::vector<Real> buildUniformRadialGrid(int nNodes, Real rMin, Real rMax);

// Build a uniform energy grid epsilon_i = E_threshold + (E_max - E_threshold)/N_E * i,
// for i = 1..N_E. Sampling starts at i=1 (never i=0) to avoid the k=0 (E=E_threshold=0)
// singularity in matchAsymptotic's matching formulas; per
// docs/planning/tise-task-breakdown.md B2, this generalizes the PDF's
// [0, E_max] grid to a user-supplied [E_threshold, E_max] window.
std::vector<Real> buildEnergyGrid(Real E_threshold, Real E_max, int N_E);

// Radial hydrogen-like potential V_L(x) = L(L+1)/(2x^2) - 1/x.
double radialPotential(double x, int L);

// True if x falls within `interval`, e.g. "[0, 20)", "(-inf, inf)". Bounds may
// be finite numbers or (+/-)inf/infinity (case-insensitive); brackets select
// inclusive ([, ]) vs. exclusive ((, )) endpoints. Throws std::runtime_error
// if `interval` doesn't match the expected syntax.
bool inInterval(double x, const std::string &interval);

// Given a piecewise potential (domain string -> muparser expression string in
// `x`), find the piece whose domain contains x and evaluate it there. Throws
// std::runtime_error if no piece's domain covers x. If more than one piece's
// domain covers x, the first match under this map's own (domain-string
// sorted, NOT declaration-order) iteration is silently used -- see
// validateNoOverlappingPotentialPieces below to guard against that
// ambiguity ahead of time, since evaluateFunction itself is called on the
// hot path (once per quadrature point) and does not re-check it.
double evaluateFunction(std::map<std::string, std::string> function, double x);

// One-time validation (call once per config/potential, NOT per
// evaluateFunction call): throws std::runtime_error naming the two
// colliding domain strings if any two pieces of `potential` cover a common
// x. Pragmatic, sampling-based check (each piece's own two boundary values
// plus its midpoint, probed against every other piece) rather than exact
// interval algebra -- sufficient for the realistic case of a finite list of
// contiguous/near-contiguous config.yaml potential pieces; does not
// separately check for GAPS (a query that lands in the gap between pieces
// is evaluateFunction's own, unrelated "does not cover x" error).
void validateNoOverlappingPotentialPieces(const std::map<std::string, std::string> &potential);

// Parse one config.yaml `potential` list entry -- a single-quoted Python
// dict-literal string, e.g. "{'domain': '(0, 100]', 'function': '-1/x'}"
// (controller.py's parse_potential_piece, via ast.literal_eval, is the
// Python-side twin of this parser -- NOT JSON, unlike main.cpp's argv-based
// parsePiecewise). Returns {domain, function}. Throws std::runtime_error on
// a missing key or malformed input; does not validate the domain/function
// strings themselves -- that happens the first time evaluateFunction uses
// them.
std::pair<std::string, std::string> parsePotentialPiece(const std::string &piece);

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
//
// === A4b: generalized drop-set (REQ-F-050 table row 4) ===
// `dropSet`, if provided, is the exact set of 1-based physical B-spline
// indices to exclude from the eigenproblem (may include interior indices,
// e.g. B-splines flanking a singular point -- see bSplinesTouchingX).
// nullopt (the default) reproduces the original hardcoded convention,
// dropping exactly {1, bs.getNBSplines()} -- bit-identical to this
// function's pre-A4b behavior (proved in docs/planning/engineer-a-plan-A4b.md).
// `nEn` MUST equal bs.getNBSplines() minus the resolved drop-set's size, or
// this throws std::runtime_error (see the plan doc for why this check
// exists -- a mismatch here previously risked a silent out-of-bounds write).
// `case3RightR`/`case3RightDelta`: Case-3 (irregular-asymptote) remediation.
// classifyAsymptote can DETECT an irregular tail and recommend a transition
// width, but detection alone doesn't change anything filled here -- passing
// both (matching `case3WindowFunction`'s R/delta, DomainSide::Right) makes
// the potential evaluate via evaluateWindowedPotential instead of plain
// evaluateFunction, smoothly tapering it to 0 over [R-delta, R] rather than
// integrating the raw (possibly irregular/divergent) tail right up to the
// wall. Left-side windowing is not supported (no current caller needs it).
// Both nullopt (the default) reproduces plain evaluateFunction, unchanged.
std::pair<std::vector<Real>, std::vector<Real>>
fillBandedMatrices(const bspline::BSpline &bs, int nEn, int order, int L,
                    std::map<std::string, std::string> potential,
                    std::optional<std::vector<int>> dropSet = std::nullopt,
                    std::optional<Real> case3RightR = std::nullopt,
                    std::optional<Real> case3RightDelta = std::nullopt);

// Given the set of BSplines, Hamiltonian, and eigenvectors, solve for:
// < phi_n | H | B_N > and < phi_n | B_N >, for each eigenvector
// These are returned as {coeffs1, coeffs2}, each with length equal to the number of eigenvectors
//
// `nBSplinesOpt`/`dropSet` generalize this beyond the classic "only B_1
// dropped" assumption (A4b, REQ-F-050): `dropSet` is the exact set of
// physical B-spline indices excluded from the bound-state basis that
// produced `Hmat`/`Smat`/`eigen` (B_N itself must NOT be in it -- its raw
// column is what this function extracts). Defaults (`nullopt`/`nullopt`)
// reproduce the original hardcoded {1}-only assumption bit-identically.
std::pair<std::vector<Real>, std::vector<Real>> precomputeBoundaryCoupling(
    int order, int nEn, std::vector<Real> Hmat, std::vector<Real> Smat, EigenResult eigen,
    std::optional<int> nBSplinesOpt = std::nullopt,
    std::optional<std::vector<int>> dropSet = std::nullopt);

// given the output of precomputeBoundaryCoupling, construct continuum states |\bar\psi_E> (linear combinations of eigenvectors)
// for each energy E on the input grid. The result is stored as a 2-D vector with (# Energy points) elements of length (# eigenvectors)
// note that each energy has a state with (# eigenvectors + 1) elements, but the coefficient for the last element is always 1 (and corresponds to B_N).
// `nBSplinesOpt`/`dropSet`: see precomputeBoundaryCoupling, forwarded as-is.
//
// The c_n(E) = (...)/(E-E_n) sum has a pole at every confined eigenvalue
// (bound OR an unfiltered box-discretization artifact, per ADR-0007) --
// `poleTolFraction`/`warnOut` add a diagnostic-only guard (no change to the
// computed values): if a grid point in `grid` lands within
// `poleTolFraction * (grid[1]-grid[0])` of any eigen.values[i] (only
// checked when grid.size() >= 2), one warning line is written to `warnOut`
// naming the energy and the closest eigenvalue. This is a physics warning
// (SDD Sec 8: computation completed, result may be unreliable), not an
// error -- the caller decides whether/how to surface it (see
// tise_solver_main.cpp's warnings.json wiring).
std::vector<std::vector<Real>> buildContinuumState(
    int order, int nEn,
    std::vector<Real> Hmat,
    std::vector<Real> Smat,
    EigenResult eigen,
    std::vector<Real> grid,
    std::optional<int> nBSplinesOpt = std::nullopt,
    std::optional<std::vector<int>> dropSet = std::nullopt,
    Real poleTolFraction = 0.1,
    std::ostream &warnOut = std::cerr
);

// `states` holds, per buildContinuumState's contract, coefficients of the
// confined eigenstates {phi_n} (plus the B_N term) -- not raw B-spline
// coefficients. `eigen` (the same EigenResult states was built from) is
// required to transform each energy's coefficients into true B-spline
// coefficients before evaluating psi_E(R) and psi_E'(R). `dropSet`: see
// precomputeBoundaryCoupling -- must match whatever drop-set actually
// produced `eigen`/`states`, or coefficients will be silently misattributed
// to the wrong physical B-spline indices before bs.eval is called.
//
// `order`/`Hmat`/`Smat`: the same inputs buildContinuumState took to build
// `states` -- needed here too, internally, to evaluate psi_E/psi_E'(R) at
// auxiliary energies a small step off the requested grid (see `fineDE`
// below), which `states` alone (built only at the production grid) can't
// provide.
//
// `result.delta` is stored as a CONTINUOUS, unwrapped trajectory (first
// point normalized into (-pi/2, pi/2], each subsequent point tracked for
// pi-periodic atan branch jumps), not raw atan(...)-kR -- see
// docs/tests/reports/8236239/free_particle.md ("phase_shifts.png ...
// uninformative as drawn"). `result.dDeltaDE[i]` is a central difference of
// this same raw-but-unwrapped delta over a small internal step `fineDE`
// around grid[i] (clamped to keep E-h > 0), NOT a finite difference across
// the (often much coarser, and non-smooth near delta=pi/4 mod pi/2 for any
// step size) production grid -- see
// docs/tests/reports/8236239/finite_square_well.md section 7 item 1 for why
// the previous sin(2*delta)/cos(2*delta) construction across the production
// grid was replaced.
AsymptoticResult matchAsymptotic(const bspline::BSpline &bs, std::vector<std::vector<Real>> states, const EigenResult &eigen, std::vector<Real> grid, Real R,
                                  int order, const std::vector<Real> &Hmat, const std::vector<Real> &Smat,
                                  std::optional<std::vector<int>> dropSet = std::nullopt,
                                  Real fineDE = 1e-3);

// Writes phase_shifts.dat-style output (epsilon_i, delta, dDeltaDE) to `out`,
// and one continuum_state_NNN.dat-style block (x, psi_E(x)) per energy to
// `stateOut[i]`. `grid` and `states` are the energy grid and per-energy
// B-spline coefficient vectors produced by buildContinuumState; `result` is
// matchAsymptotic's output for that same grid.
//
// `states` holds, per buildContinuumState's contract, coefficients of the
// confined eigenstates {phi_n} (plus the B_N term) -- not raw B-spline
// coefficients (see matchAsymptotic's doc comment above, which this
// function's `eigen`/`dropSet` parameters exist to satisfy identically):
// `eigen` (the same EigenResult `states` was built from) is required to
// transform each energy's coefficients into true B-spline coefficients
// before evaluating psi_E(x). `dropSet` must match whatever drop-set
// actually produced `eigen`/`states`, or coefficients will be silently
// misattributed to the wrong physical B-spline indices before bs.eval is
// called.
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
                     std::optional<std::vector<int>> dropSet = std::nullopt);

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
// Singular/Continuous produce no knot -- Singular's remediation is B-spline
// removal (A4b, `bSplinesTouchingX` below), a different function from this
// one's knot-insertion job, not an unimplemented gap: A4b is implemented
// and, since ADR-0009, wired into both solveTISE and tise_solver_main.cpp
// via buildStrategicGridAndDropSet.
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
// them clamp to the same boundary knot), not just one -- see
// docs/planning/engineer-a-plan-A4b.md for the worked example.
std::vector<int> bSplinesTouchingX(int nNodes, int order, const std::vector<Real> &grid, Real x);

// Analytic hydrogenic energy: E = -1 / (2 * (n + L)^2).
Real analyticHydrogenEnergy(int n, int L);

// Difference between computed eigenvalue and analytic energy for state n.
Real eigenvalueError(Real computed, int n, int L);

// Extract eigenvector column iEn (1-based) from the column-major evec array
// and embed it into a zero-padded vector of length nBSplines. Physical
// indices in `dropSet` (or {1,nBSplines} if omitted -- see
// fillBandedMatrices, same contract) get 0; kept indices get their
// eigenvector component, in ascending physical-index order (must match
// whatever dropSet fillBandedMatrices was called with, or coefficients
// will be silently misattributed -- both functions share the same
// resolveDropSet helper, so passing the same dropSet to both guarantees
// consistency).
std::vector<Real> eigenstateCoefficients(const std::vector<Real> &evec,
                                          int iEn,
                                          int nEn,
                                          int nBSplines,
                                          std::optional<std::vector<int>> dropSet = std::nullopt);

// Write a single eigenstate to `out`: npts lines of "x  psi(x)".
void writeEigenstate(std::ostream &out,
                     const bspline::BSpline &bs,
                     const std::vector<Real> &coeffs,
                     int npts,
                     Real rMin,
                     Real rMax);

// Write eigenvalues.dat: 0-based index, E_n, one line per state, for the
// first nStates entries of er.values (ascending, per EigenResult's own
// contract). Per ADR-0007, no bound/continuum filtering is applied here --
// callers pass er.dim to write every computed state.
void writeEigenvalues(std::ostream &out, const EigenResult &er, int nStates);

// Write eigenvectors.dat: nBSplines rows x nStates columns. Column j is the
// full, zero-padded B-spline coefficient vector for eigenstate j+1
// (1-based, per eigenstateCoefficients). `dropSet` MUST match whatever
// drop-set er was actually diagonalized under (defaults to the classic
// {1,nBSplines} convention, matching fillBandedMatrices(..., nEn+1, ...,
// {1}) truncated to its leading nEn columns for the solve -- see
// docs/SDD.md §6.3) -- passing a mismatched dropSet silently misattributes
// coefficients, same contract as eigenstateCoefficients itself.
void writeEigenvectors(std::ostream &out, const EigenResult &er, int nBSplines, int nStates,
                        std::optional<std::vector<int>> dropSet = std::nullopt);

// Write hamiltonian.dat/overlap.dat: preserves fillBandedMatrices' own
// column-major banded layout (order rows x nEn cols, element (row,col) at
// mat[(row-1)+(col-1)*order]) as plain text, row-major. Reads only the
// leading order*nEn elements of `mat` -- a caller may pass a wider array
// (e.g. the order*(nEn+1) continuum-coupling fill) and get back just the
// nEn x nEn block that was actually diagonalized.
void writeBandedMatrix(std::ostream &out, const std::vector<Real> &mat,
                        int order, int nEn, const std::string &description);

// Bundles solveTISE's full output: the diagonalized EigenResult plus every
// piece of basis-construction context a caller needs to correctly decode it
// downstream (docs/planning/engineer-a-plan-A4-wiring.md, Gap 3). Before
// this, solveTISE returned only EigenResult, forcing callers (main.cpp) to
// independently reconstruct the grid/BSpline/dropSet -- harmless only while
// solveTISE always built a plain uniform grid with the classic
// {1,nBSplines} drop-set. Now that solveTISE may build a strategic
// (non-uniform) grid and/or a non-classic drop-set (REQ-F-050), an
// independent reconstruction would silently diverge from what was actually
// solved.
// === ADR-0009: shared grid+drop-set construction (unifies solveTISE/tise_solver) ===
// Extracted from solveTISE's own grid/dropset construction (was previously
// inlined there only) so tise_solver_main.cpp -- which used to build a
// plain uniform grid and hardcode the classic {1} drop-set, structurally
// unable to reach strategic node placement (REQ-F-050) or A4b interior
// singular-B-spline removal -- can share the exact same logic instead of a
// second, divergent copy (the risk ADR-0008 itself named as its revisit
// trigger: "if the two copies drift and cause a bug"). Automatically
// strategic if the potential has detectable Step/StitchedKink/Singular
// structure; a potential with none produces an unchanged uniform grid and
// the classic {1} drop-set, byte-identical to the pre-ADR-0009 behavior.
struct StrategicGridResult
{
    std::vector<Real> grid;    // the exact (possibly non-uniform/strategic) physical
                                // grid passed to bs.init
    bspline::BSpline bs;       // the exact basis constructed from `grid`
    int nBSplines;              // == bs.getNBSplines(); duplicated for convenience
    std::vector<int> fillDropSet; // for fillBandedMatrices -- B_N (nBSplines) is
                                // deliberately never included (its raw H/S column is
                                // needed by continuum construction, not dropped from
                                // the fill; it's excluded from the DIAGONALIZATION by
                                // the nEnFilled/nEnBound truncation below instead)
    int nEnBound;                // == nBSplines - fillDropSet.size() - 1
    bool rightEdgeSingular;    // true iff the potential is singular at x=rMax --
                                // matchAsymptotic's flat-asymptote assumption there
                                // is affected regardless of B-spline removal
};

StrategicGridResult buildStrategicGridAndDropSet(int nNodes, int order, Real rMin, Real rMax,
                                                   const std::map<std::string, std::string> &potential);

struct SolveTISEResult
{
    EigenResult eigen;         // as returned by solveGeneralizedEigenproblem
    bspline::BSpline bs;       // the exact basis diagonalized against -- reuse
                                // directly for eigenstateCoefficients/writeEigenstate/
                                // runTimeEvolution; do not rebuild
    std::vector<Real> grid;    // the exact (possibly non-uniform/strategic) physical
                                // grid passed to bs.init -- needed by minInterNodeGap
    int nBSplines;              // == bs.getNBSplines(); duplicated for convenience
    std::vector<int> dropSet;  // FULL set of physical B-spline indices excluded from
                                // the bound-state basis (classic wall unioned with any
                                // A4b Singular-join clusters). Pass directly as
                                // eigenstateCoefficients' dropSet argument -- eigen.vectors
                                // was diagonalized against exactly this exclusion set.
};

// Top-level TISE solver: build grid (automatically strategic per REQ-F-050
// if the potential has detectable Step/StitchedKink/Singular structure --
// see docs/planning/engineer-a-plan-A4-wiring-design.md), fill matrices
// (with singular-join B-splines removed per A4b), diagonalise, then
// construct and write continuum states/phase shifts on the energy grid
// [E_threshold, E_max] (N_E points, per buildEnergyGrid).
// `potential` is the piecewise potential passed through to fillBandedMatrices.
SolveTISEResult solveTISE(int nNodes, int order, Real rMin, Real rMax, int L, std::map<std::string, std::string> potential,
                           Real E_threshold, Real E_max, int N_E);

// === A5: E_acc continuum-accuracy warning (REQ-F-040, warning half) ===
// Reduce a (possibly non-uniform, possibly containing degenerate/repeated
// knots from buildStrategicRadialGrid) grid to the minimum spacing between
// DISTINCT physical node locations, for use as computeEAcc's nodeSpacing on
// a strategic grid. Repeated copies of the same x (the degenerate knots
// strategic placement inserts at Step/StitchedKink joins) do not represent
// a second, vanishingly-close physical point -- they increase local knot
// multiplicity at one location -- so they're collapsed before taking the
// minimum gap; naively taking min(grid[i+1]-grid[i]) over a raw strategic
// grid would otherwise return 0.0 at any such knot. `grid` must be sorted
// non-decreasing. Throws std::runtime_error if fewer than 2 distinct points
// remain after collapsing.
Real minInterNodeGap(const std::vector<Real> &grid, Real tol = 1e-12);

// Basis accuracy ceiling (docs/planning/architecture-06-20.md "Continuum
// range" -> "Basis accuracy limit"; duplicated docs/SDD.md Appendix B
// "Continuum range" entry; see also SDD Sec. 6.4, Sec. 8). States whose
// half de Broglie wavelength is commensurate with or smaller than the
// B-spline node spacing cannot be accurately represented:
//   lambda/2 = pi/k <~ dx_node  ==>  k >~ pi/dx_node  ==>  E >~ pi^2 / (2 m dx_node^2)
// The source relation is an asymptotic ("<~"/">~") scaling bound, not an
// exact equality; this function takes its leading-order coefficient as
// the concrete E_acc threshold, per this task's closed-form "Done when"
// criterion -- treat the result as an order-of-magnitude ceiling, not a
// razor-sharp cutoff. `nodeSpacing` is a single scalar spacing; for a
// non-uniform grid (see A4's buildStrategicRadialGrid), pass
// minInterNodeGap(grid) above, not an average.
Real computeEAcc(Real nodeSpacing, Real mass);

// Warns (to warnOut, default stderr -- SDD Sec. 8's "physics warning"
// class: computation completed, but a result may be unreliable) iff the
// requested continuum ceiling eMax exceeds the basis accuracy ceiling
// eAcc (from computeEAcc), and returns true iff it did so. Strict '>':
// eMax exactly equal to eAcc is the marginal case and does not warn,
// mirroring checkWellContainment's/classifyBoundStates' existing
// strict-inequality boundary convention. Standalone: no dependency on
// Engineer B's energy-grid loop (B2) -- wiring this at that call site is
// future integration work, matching A1-A4's "standalone, unwired"
// precedent.
bool warnIfContinuumExceedsEAcc(Real eMax, Real eAcc, std::ostream &warnOut = std::cerr);

} // namespace tise
