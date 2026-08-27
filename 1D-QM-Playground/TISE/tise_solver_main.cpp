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

        // B-spline basis parameters.
        const int nNodes = config["bspline"]["n_nodes"].as<int>();
        const int order  = config["bspline"]["order"].as<int>();
        const tise::Real rMin = config["bspline"]["domain"][0].as<tise::Real>();
        const tise::Real rMax = config["bspline"]["domain"][1].as<tise::Real>();

        // Piecewise potential.
        const std::map<std::string, std::string> potential = parsePotentialConfig(config["potential"]);

        // Construct the B-spline basis -- mirrors tise::solveTISE's own construction.
        const std::vector<tise::Real> grid = tise::buildUniformRadialGrid(nNodes, rMin, rMax);
        bspline::BSpline bs;
        if (bs.init(nNodes, order, grid) != 0)
            throw std::runtime_error("BSpline::init failed");

        // Basis size / eigenproblem dimension -- matches solveTISE's own pattern.
        const int nBSplines = bs.getNBSplines();
        const int nEn       = nBSplines - 2;

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

        // Fill banded matrices. dropSet={1}, nEn+1 columns: continuum-ready
        // (column nEn+1 is B_N's column, needed by the continuum path
        // below). Truncating to the leading nEn columns for the solve
        // (next step) reproduces the classic {1,nBSplines} default drop-set
        // exactly -- verified bit-identical during this feature's planning,
        // not an approximation. L is fixed at 0: fillBandedMatrices no
        // longer uses it to select the potential (only eigenvalueError's
        // irrelevant hydrogen-analytic comparison does), and config.yaml
        // has no L field -- centrifugal terms are baked into the potential
        // expression itself.
        auto [H, S] = tise::fillBandedMatrices(bs, nEn + 1, order, /*L=*/0, potential, std::vector<int>{1});

        // Solve. H, S are passed by value -- solveGeneralizedEigenproblem's
        // internal LAPACK call overwrites its own copies, not these, so H/S
        // remain valid below for hamiltonian.dat/overlap.dat and the
        // continuum path.
        tise::EigenResult er = tise::solveGeneralizedEigenproblem(H, S, nEn, order);

        std::vector<WarningEntry> warnings;

        // Case-3-asymptote warning: classifyAsymptote's own documented
        // precondition is "the caller already knows this side is
        // unbounded." Approximate that by checking whether the potential's
        // own domain coverage actually extends far beyond the box -- for a
        // potential capped exactly at the box edge (the common case, e.g.
        // the real config.yaml), there is no unbounded side to classify and
        // this is skipped entirely.
        {
            const tise::Real probeX = rMax + 1.0e6 * std::max(std::abs(rMax), tise::Real(1.0));
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
                    tise::classifyAsymptote(potential, tise::SpatialDomain{rMin, rMax},
                                             tise::DomainSide::Right, warnOut);
                    if (!warnOut.str().empty())
                        addWarning(warnings, "physics", warnOut.str());
                }
                catch (const std::exception &)
                {
                    // Classification not applicable for this potential/domain
                    // shape -- treated as "no diagnostic available," not a
                    // solver failure.
                }
            }
        }

        // Continuum construction, gated on tise.continuum.enabled.
        if (continuumEnabled)
        {
            // mass=1.0 hardcoded, matching fillBandedMatrices' own
            // internal kinetic-energy term (which already hardcodes /2.0,
            // i.e. mass=1 baked into the matrix fill itself) -- reading
            // config["physics"]["mass"] only here would suggest it's
            // configurable when the core solve ignores it entirely.
            const tise::Real nodeSpacing = (rMax - rMin) / static_cast<tise::Real>(nNodes - 1);
            const tise::Real eAcc = tise::computeEAcc(nodeSpacing, /*mass=*/1.0);
            {
                std::ostringstream warnOut;
                if (tise::warnIfContinuumExceedsEAcc(E_max, eAcc, warnOut))
                    addWarning(warnings, "physics", warnOut.str());
            }

            auto energyGrid = tise::buildEnergyGrid(E_threshold, E_max, n_energies);
            auto states = tise::buildContinuumState(order, nEn, H, S, er, energyGrid);
            auto ar = tise::matchAsymptotic(bs, states, er, energyGrid, /*R=*/rMax);

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
                                      n_pts, rMin, rMax);
        }

        // Core output files -- all nEn states, per ADR-0007.
        {
            std::ofstream out(outputDir / "eigenvalues.dat");
            tise::writeEigenvalues(out, er, er.dim);
        }
        {
            std::ofstream out(outputDir / "eigenvectors.dat");
            tise::writeEigenvectors(out, er, nBSplines, er.dim);
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
