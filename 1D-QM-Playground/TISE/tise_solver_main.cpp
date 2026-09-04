// tise_solver_main.cpp
//
// CLI for the `tise_solver` executable described in docs/SDD.md §7.2.1
// (Controller to TISE Solver). Parses config.yaml's bspline/potential/tise
// blocks, solves the bound-state (and, if enabled, continuum) problem via
// the real tise:: library, and writes the data/tise/ output contract
// (docs/SDD.md §6.3). Per ADR-0007, eigenvalues.dat/eigenvectors.dat
// contain all computed states -- no bound-state filtering is applied here.
//
// Usage:
//   tise_solver --config <config.yaml> --output-dir <output_dir>

#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>
#include <yaml-cpp/yaml.h>

#include "BSpline.hpp"
#include "tise.hpp"

namespace
{

constexpr const char *USAGE = "usage: tise_solver --config <path> --output-dir <path>";

struct CliArgs
{
    std::string configPath;
    std::string outputDir;
};

// Parse the two required flags. Anything unexpected -- missing flags,
// unrecognized arguments, or a flag with no following value -- is a usage
// error and throws std::runtime_error.
CliArgs parseArgs(int argc, char *argv[])
{
    std::string configPath;
    std::string outputDir;
    bool haveConfig = false;
    bool haveOutputDir = false;

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];

        if (arg == "--config")
        {
            if (i + 1 >= argc)
                throw std::runtime_error(std::string("--config requires a value\n") + USAGE);
            configPath = argv[++i];
            haveConfig = true;
        }
        else if (arg == "--output-dir")
        {
            if (i + 1 >= argc)
                throw std::runtime_error(std::string("--output-dir requires a value\n") + USAGE);
            outputDir = argv[++i];
            haveOutputDir = true;
        }
        else
        {
            throw std::runtime_error("unrecognized argument '" + arg + "'\n" + USAGE);
        }
    }

    if (!haveConfig || !haveOutputDir)
        throw std::runtime_error(std::string("missing required argument(s)\n") + USAGE);

    return CliArgs{configPath, outputDir};
}

// Validate and extract config["tise"]["continuum"]["n_energies"]. Deliberately
// separated from the continuum-output path and called BEFORE any output file
// is written in main(): the Controller<->TISE contract (docs/SDD.md §7.2.1)
// treats any non-zero exit as a "hard failure" that must leave no partial
// output behind, so everything that could invalidate this run has to be
// validated before the first byte of any output file is written.
int validateNEnergies(const YAML::Node &continuumNode)
{
    // Missing or wrong-typed n_energies throws YAML::TypedBadConversion<int>
    // (or YAML::InvalidNode if the key is absent entirely) -- both derive
    // from std::exception and propagate to main()'s catch block unmodified.
    const int nEnergies = continuumNode["n_energies"].as<int>();
    if (nEnergies <= 0)
        throw std::runtime_error("tise.continuum.n_energies must be a positive integer");
    return nEnergies;
}

// Loop over config["potential"]'s YAML sequence, parsing each entry via
// tise::parsePotentialPiece (the yaml-cpp-independent half of this parser,
// unit-tested in TISE/tests/test_tise.cpp). This thin wrapper is the only
// yaml-cpp-touching potential-parsing code, kept local to this file so
// tise.hpp stays free of the yaml-cpp coupling tise_solver's CMake scoping
// deliberately avoids (H-BoundStates must not need yaml-cpp).
std::map<std::string, std::string> parsePotentialConfig(const YAML::Node &potentialNode)
{
    if (!potentialNode || !potentialNode.IsSequence() || potentialNode.size() == 0)
        throw std::runtime_error("'potential' must be a non-empty YAML list of piece strings");

    std::map<std::string, std::string> potential;
    for (const auto &pieceNode : potentialNode)
    {
        auto [domain, function] = tise::parsePotentialPiece(pieceNode.as<std::string>());
        potential[domain] = function;
    }
    return potential;
}

struct WarningEntry
{
    std::string category;
    std::string message;
};

// Records one warning: SDD §8 requires warnings reach both stderr and the
// warnings.json sidecar.
void addWarning(std::vector<WarningEntry> &warnings, const std::string &category, std::string text)
{
    while (!text.empty() && std::isspace(static_cast<unsigned char>(text.back())))
        text.pop_back();
    std::cerr << "tise_solver: warning: " << text << "\n";
    warnings.push_back({category, std::move(text)});
}

} // namespace

int main(int argc, char *argv[])
{
    try
    {
        const CliArgs args = parseArgs(argc, argv);

        // YAML::LoadFile throws YAML::BadFile if configPath doesn't exist,
        // or YAML::ParserException if it's malformed. Let both propagate to
        // the catch block below rather than handling them here.
        const YAML::Node config = YAML::LoadFile(args.configPath);

        std::filesystem::create_directories(args.outputDir);
        const std::filesystem::path outputDir(args.outputDir);

        // physics.mass/physics.hbar guard-rail: both fields are documented
        // in config.yaml but consumed nowhere -- fillBandedMatrices
        // hardcodes mass=1 internally (matching computeEAcc's own
        // mass=1.0 call and k=sqrt(2*E) throughout), so a user setting
        // either to anything else previously got silent wrong physics with
        // no error. Guard-rail only, not full generalization (that would
        // touch the kinetic-energy term/k=sqrt(2E)/computeEAcc throughout)
        // -- an honest config error instead.
        if (config["physics"])
        {
            if (config["physics"]["mass"] && config["physics"]["mass"].as<tise::Real>() != 1.0)
                throw std::runtime_error("physics.mass is fixed at 1.0 internally (fillBandedMatrices' "
                                          "kinetic-energy term/computeEAcc/k=sqrt(2E) all hardcode it); "
                                          "setting it to anything else is not yet supported -- remove the "
                                          "field or set it to 1.0.");
            if (config["physics"]["hbar"] && config["physics"]["hbar"].as<tise::Real>() != 1.0)
                throw std::runtime_error("physics.hbar is fixed at 1.0 internally (atomic units throughout "
                                          "this solver); setting it to anything else is not yet supported -- "
                                          "remove the field or set it to 1.0.");
        }

        // B-spline basis parameters.
        const int nNodes = config["bspline"]["n_nodes"].as<int>();
        const int order  = config["bspline"]["order"].as<int>();
        const tise::Real rMin = config["bspline"]["domain"][0].as<tise::Real>();
        const tise::Real rMax = config["bspline"]["domain"][1].as<tise::Real>();

        // Piecewise potential.
        const std::map<std::string, std::string> potential = parsePotentialConfig(config["potential"]);

        // One-time overlap validation (not per-evaluateFunction-call): a
        // malformed config with two pieces covering a common x previously
        // resolved silently via std::map's own domain-string sort order
        // (not declaration order) with no error at all -- now throws a
        // clear, named error instead. controller.py's own
        // validate_potential_tiling already checks this on the Python
        // side, but tise_solver is independently invocable (SDD's own
        // architecture) and previously had no such guard of its own.
        tise::validateNoOverlappingPotentialPieces(potential);

        // Construct the B-spline basis: automatically strategic per REQ-F-050
        // if the potential has detectable Step/StitchedKink/Singular
        // structure, with A4b interior-singular-B-spline removal applied --
        // shared with H-BoundStates'/solveTISE's own construction (ADR-0009;
        // this file previously built a plain uniform grid + hardcoded
        // classic {1} drop-set only, structurally unable to reach either --
        // see docs/planning/tise-release-readiness-plan.md Part A). A
        // potential with no detectable structure produces the same uniform
        // grid + {1} drop-set as before, byte-identical.
        tise::StrategicGridResult sgr = tise::buildStrategicGridAndDropSet(nNodes, order, rMin, rMax, potential);
        bspline::BSpline &bs = sgr.bs;
        const int nBSplines = sgr.nBSplines;
        const int nEn       = sgr.nEnBound;

        // Spatial grid density for eigenstate_NNN.dat, written unconditionally
        // below regardless of tise.continuum.enabled -- previously read by
        // nothing in this file.
        const int nPtsEigenstate = config["tise"]["n_pts_eigenstate"].as<int>();

        // Continuum settings. yaml-cpp idiom: an undefined/missing Node
        // converts to false without throwing, so configs that omit
        // tise.continuum entirely, or set enabled: false, correctly skip
        // the branches below.
        const bool continuumEnabled = config["tise"] && config["tise"]["continuum"] &&
                                       config["tise"]["continuum"]["enabled"] &&
                                       config["tise"]["continuum"]["enabled"].as<bool>();

        // Validate everything this run will need BEFORE writing anything
        // (see validateNEnergies' comment: a non-zero exit must leave no
        // partial output on disk).
        tise::Real E_threshold = 0.0, E_max = 0.0;
        int n_energies = 0, n_pts = 0;
        if (continuumEnabled)
        {
            const YAML::Node continuumNode = config["tise"]["continuum"];
            n_energies  = validateNEnergies(continuumNode);
            E_threshold = continuumNode["E_threshold"].as<tise::Real>();
            E_max       = continuumNode["E_max"].as<tise::Real>();
            n_pts       = continuumNode["n_pts"].as<int>();
        }

        std::vector<WarningEntry> warnings;

        // Case-3-asymptote classification, BEFORE fillBandedMatrices (not
        // after, as originally): an Irregular classification's
        // recommendedTransitionWidth needs to reach fillBandedMatrices so
        // it can actually taper the potential there (previously this
        // classification only produced a warning -- detection without
        // remediation, see docs/planning/tise-release-readiness-plan.md
        // Part B). classifyAsymptote's own documented precondition is "the
        // caller already knows this side is unbounded." Approximate that
        // by checking whether the potential's own domain coverage actually
        // extends far beyond the box -- for a potential capped exactly at
        // the box edge (the common case, e.g. the real config.yaml), there
        // is no unbounded side to classify and this is skipped entirely.
        std::optional<tise::Real> case3RightR, case3RightDelta;
        {
            // Not a shared/parameterized library function -- there is only
            // this one caller, so these stay plain local constants rather
            // than optional parameters on some shared helper (see the
            // comment on kUnboundedProbeMagnitude in
            // validateNoOverlappingPotentialPieces for why this file, that
            // function, and classifyAsymptote's own internal `scale` each
            // use their own independent "probe far away" constant instead
            // of a common one -- each solves a distinct problem).
            constexpr tise::Real kProbeDistanceMultiplier = 1.0e6; // how many
                // box-widths beyond rMax to probe -- far enough that a
                // potential piece genuinely extending to +inf is virtually
                // certain to still cover this point, without risking
                // overflow the way an actual +inf-magnitude probe would.
            constexpr tise::Real kProbeDistanceFloor = 1.0; // floor on
                // std::abs(rMax) so a box edge at/near rMax=0 still gets a
                // meaningfully far-away probe point instead of one at x=0.
            const tise::Real probeX = rMax + kProbeDistanceMultiplier *
                                       std::max(std::abs(rMax), kProbeDistanceFloor);
            bool coveredBeyondDomain = false;
            for (const auto &[domainStr, fn] : potential)
            {
                (void)fn;
                if (tise::inInterval(probeX, domainStr))
                {
                    coveredBeyondDomain = true;
                    break;
                }
            }
            if (coveredBeyondDomain)
            {
                try
                {
                    std::ostringstream warnOut;
                    auto classification = tise::classifyAsymptote(
                        potential, tise::SpatialDomain{rMin, rMax}, tise::DomainSide::Right, warnOut);
                    if (!warnOut.str().empty())
                        addWarning(warnings, "physics", warnOut.str());
                    if (classification.asymptoteCase == tise::AsymptoteCase::Irregular)
                    {
                        std::ostringstream msg;
                        msg << "potential's right-edge tail is Irregular (Case 3)";
                        // Gated on continuum.enabled (release-readiness
                        // follow-up plan, Part C): the taper exists solely
                        // to make matchAsymptotic's flat-asymptote
                        // assumption valid at R=rMax for continuum
                        // matching. Applying it unconditionally shifted
                        // bound-state energies by up to 1.5e-4 Ha for zero
                        // benefit whenever continuum was disabled
                        // (docs/tests/reports/8236239/case3_irregular_tail.md).
                        if (continuumEnabled)
                        {
                            case3RightR = rMax;
                            case3RightDelta = classification.recommendedTransitionWidth;
                            msg << "; tapering it to 0 over a transition width of "
                                << classification.recommendedTransitionWidth
                                << " before x=" << rMax << " rather than integrating the raw tail up "
                                   "to the wall (see evaluateWindowedPotential).";
                        }
                        else
                        {
                            msg << "; continuum is disabled, so the raw (untapered) potential is "
                                   "integrated to the wall -- the bound-state spectrum is unaffected.";
                        }
                        addWarning(warnings, "physics", msg.str());
                    }
                }
                catch (const std::exception &)
                {
                    // Classification not applicable for this potential/domain
                    // shape -- treated as "no diagnostic available," not a
                    // solver failure.
                }
            }
        }

        // Fill banded matrices. sgr.fillDropSet, nEn+1 columns: continuum-ready
        // (column nEn+1 is B_N's column, needed by the continuum path
        // below). Truncating to the leading nEn columns for the solve
        // (next step) reproduces sgr.fillDropSet's exclusion exactly. L is
        // fixed at 0: fillBandedMatrices no longer uses it to select the
        // potential (only eigenvalueError's irrelevant hydrogen-analytic
        // comparison does), and config.yaml has no L field -- centrifugal
        // terms are baked into the potential expression itself. case3RightR/
        // Delta (set above, nullopt for the common case) taper the
        // potential near the wall when Case 3 was detected.
        auto [H, S] = tise::fillBandedMatrices(bs, nEn + 1, order, /*L=*/0, potential, sgr.fillDropSet,
                                                case3RightR, case3RightDelta);

        // Solve. H, S are passed by value -- solveGeneralizedEigenproblem's
        // internal LAPACK call overwrites its own copies, not these, so H/S
        // remain valid below for hamiltonian.dat/overlap.dat and the
        // continuum path.
        tise::EigenResult er = tise::solveGeneralizedEigenproblem(H, S, nEn, order);

        // classifyBoundStates (REQ-F-020): informational count of how many
        // of the nEn computed states are below the ionization threshold
        // E=0.0. Purely informational, not a filter -- eigenvalues.dat/
        // eigenvectors.dat still contain every computed state unfiltered,
        // per ADR-0007; this just finally surfaces the split that filter
        // deliberately doesn't apply, previously invisible anywhere in this
        // binary's output.
        {
            auto classification = tise::classifyBoundStates(er, /*threshold=*/0.0);
            std::ostringstream msg;
            msg << classification.nBound << " of " << er.dim << " computed states are "
                   "below E=0.0 (informational bound-state count; eigenvalues.dat/"
                   "eigenvectors.dat contain all computed states unfiltered, per ADR-0007).";
            addWarning(warnings, "physics", msg.str());
        }

        // Distinct from the Case-3 check above: sgr.rightEdgeSingular fires
        // when a potential piece is itself singular exactly AT x=rMax (e.g.
        // "1/(rMax-x)"), regardless of whether anything is defined beyond
        // the box -- previously only solveTISE (the H-BoundStates path)
        // could detect this, via its own bare-cerr warning; this file had
        // no access to it at all before ADR-0009's shared construction.
        // Only meaningful when continuum construction is actually going to
        // run below -- the bound-state solve itself is entirely unaffected
        // by a right-edge singularity (the same classic B_N wall-exclusion
        // that already handles hydrogen's origin singularity applies here
        // too), so there is nothing to warn about when continuum is off.
        if (sgr.rightEdgeSingular && continuumEnabled)
            addWarning(warnings, "physics",
                       "potential is singular at the right domain edge x=" + std::to_string(rMax) +
                       "; continuum phase-shift matching (matchAsymptotic) assumes a regular "
                       "boundary there, which does not hold here -- skipping continuum "
                       "construction entirely (no phase_shifts.dat/continuum_state_NNN.dat "
                       "written). The bound-state solve above is unaffected and still valid.");

        // Continuum construction, gated on tise.continuum.enabled AND not
        // singular at the right domain edge (known-solution-verification
        // follow-up plan, Part D): a right-edge singularity requires
        // psi_E(rMax)=0 physically, but B_N -- the "escape" function
        // continuum construction deliberately keeps -- is exactly the one
        // basis function that's non-zero there. The resulting "continuum
        // state" is essentially B_N alone, not a solution of the problem
        // (docs/tests/reports/8236239/right_edge_singularity.md section
        // 4.3: 30x larger at the wall than anywhere in the interior).
        // Refuse rather than write known-wrong output that a downstream
        // plot would happily draw.
        if (continuumEnabled && !sgr.rightEdgeSingular)
        {
            // mass=1.0 hardcoded, matching fillBandedMatrices' own
            // internal kinetic-energy term (which already hardcodes /2.0,
            // i.e. mass=1 baked into the matrix fill itself) -- reading
            // config["physics"]["mass"] only here would suggest it's
            // configurable when the core solve ignores it entirely.
            // minInterNodeGap (not a flat (rMax-rMin)/(nNodes-1) average):
            // sgr.grid may now be a non-uniform strategic grid, and the
            // physically-correct nodeSpacing for a non-uniform grid is its
            // minimum inter-node gap, not an average across it -- see
            // minInterNodeGap's own doc comment.
            const tise::Real nodeSpacing = tise::minInterNodeGap(sgr.grid);
            const tise::Real eAcc = tise::computeEAcc(nodeSpacing, /*mass=*/1.0);
            {
                std::ostringstream warnOut;
                if (tise::warnIfContinuumExceedsEAcc(E_max, eAcc, warnOut))
                    addWarning(warnings, "physics", warnOut.str());
            }

            auto energyGrid = tise::buildEnergyGrid(E_threshold, E_max, n_energies);
            std::ostringstream poleWarnOut;
            auto states = tise::buildContinuumState(order, nEn, H, S, er, energyGrid,
                                                      sgr.nBSplines, sgr.fillDropSet, 0.1, poleWarnOut);
            if (!poleWarnOut.str().empty())
                addWarning(warnings, "physics", poleWarnOut.str());
            auto ar = tise::matchAsymptotic(bs, states, er, energyGrid, /*R=*/rMax, order, H, S, sgr.fillDropSet);

            std::ofstream phaseShiftsOut(outputDir / "phase_shifts.dat");
            std::vector<std::ofstream> continuumStateFiles;
            continuumStateFiles.reserve(energyGrid.size());
            std::vector<std::ostream *> continuumStateOut;
            for (std::size_t i = 0; i < energyGrid.size(); ++i)
            {
                std::ostringstream oss;
                oss << "continuum_state_" << std::setw(3) << std::setfill('0') << (i + 1) << ".dat";
                continuumStateFiles.emplace_back(outputDir / oss.str());
                continuumStateOut.push_back(&continuumStateFiles.back());
            }
            tise::writeContinuumInfo(phaseShiftsOut, bs, ar, energyGrid, states, continuumStateOut,
                                      n_pts, rMin, rMax, er, sgr.fillDropSet);
        }

        // er.vectors are nEn(=nEnBound)-dimensional -- excluded from them is
        // sgr.fillDropSet's physical indices PLUS B_N (nBSplines), which was
        // truncated away separately by solveGeneralizedEigenproblem(H, S,
        // nEnBound, order) above, not by fillDropSet itself (fillBandedMatrices
        // was filled with nEn+1 columns precisely to keep B_N's column for
        // the continuum path). eigenstateCoefficients/writeEigenvectors zero-pad
        // a BOUND eigenvector (no B_N "+1" convention, unlike
        // buildContinuumState/matchAsymptotic/writeContinuumInfo above) --
        // they need this FULL exclusion set, or they'd expect one more kept
        // column than er.vectors actually has and silently misattribute
        // every coefficient from B_N's physical index onward. Mirrors
        // solveTISE's own identical fillDropSet+nBSplines construction of
        // SolveTISEResult.dropSet.
        std::vector<int> fullDropSet = sgr.fillDropSet;
        fullDropSet.push_back(sgr.nBSplines);
        std::sort(fullDropSet.begin(), fullDropSet.end());

        // Core output files -- all nEn states, per ADR-0007.
        {
            std::ofstream out(outputDir / "eigenvalues.dat");
            tise::writeEigenvalues(out, er, er.dim);
        }
        {
            std::ofstream out(outputDir / "eigenvectors.dat");
            tise::writeEigenvectors(out, er, nBSplines, er.dim, fullDropSet);
        }
        // eigenstate_NNN.dat: spatial phi_n(x), one file per bound state,
        // unconditionally (same "TISE writes everything" ADR-0007 policy as
        // eigenvalues.dat/eigenvectors.dat above -- gating what to PLOT from
        // this is analysis.py's visualization.eigenstates job, not this
        // binary's). Mirrors H-BoundStates' own main.cpp:114-123 pattern,
        // minus that demo's hydrogen-specific accuracy-threshold early exit.
        for (int j = 1; j <= nEn; ++j)
        {
            std::ostringstream oss;
            oss << "eigenstate_" << std::setw(3) << std::setfill('0') << j << ".dat";
            std::ofstream out(outputDir / oss.str());
            auto coeffs = tise::eigenstateCoefficients(er.vectors, j, nEn, nBSplines, fullDropSet);
            tise::writeEigenstate(out, bs, coeffs, nPtsEigenstate, rMin, rMax);

            // Well-containment diagnostic (A3, SDD Sec 5.2.3/6.4): a state
            // confined within the box should decay to numerically-zero
            // slope at the outer wall; a nonzero psi'(rMax) means it's
            // colliding with the wall and its energy/wavefunction may be
            // inaccurate due to box truncation -- previously computed
            // nowhere in this binary despite being fully implemented and
            // tested. Restricted to states classifyBoundStates above
            // actually calls bound (E<0.0): an unbound/continuum-like
            // state is EXPECTED to have non-negligible amplitude/slope at
            // the wall (it's not supposed to be exponentially localized in
            // the first place), so applying this check there would just
            // flag every such state unconditionally -- noise, not signal.
            if (er.values[j - 1] < 0.0)
            {
                auto containment = tise::checkWellContainment(bs, coeffs, rMax);
                if (containment.notWellContained)
                {
                    std::ostringstream msg;
                    msg << "bound state " << j << " (E=" << er.values[j - 1] << ") appears to be "
                           "colliding with the outer wall at x=" << rMax << " (psi'(rMax)="
                        << containment.psiPrimeAtBoundary << ", exceeds tolerance) -- its "
                           "energy/wavefunction may be inaccurate due to box truncation.";
                    addWarning(warnings, "physics", msg.str());
                }
            }
        }
        {
            std::ofstream out(outputDir / "hamiltonian.dat");
            tise::writeBandedMatrix(out, H, order, nEn, "hamiltonian.dat: H matrix (banded)");
        }
        {
            std::ofstream out(outputDir / "overlap.dat");
            tise::writeBandedMatrix(out, S, order, nEn, "overlap.dat: S matrix (banded)");
        }
        {
            nlohmann::json warningsJson = nlohmann::json::array();
            for (const auto &w : warnings)
                warningsJson.push_back({{"category", w.category}, {"message", w.message}});
            std::ofstream(outputDir / "warnings.json") << warningsJson.dump(2) << "\n";
        }

        // No tevol::runTimeEvolution call anywhere in this file -- TDSE
        // stays entirely out of tise_solver.
    }
    catch (const std::exception &e)
    {
        std::cerr << "tise_solver: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
