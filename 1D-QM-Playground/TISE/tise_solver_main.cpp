// tise_solver_main.cpp
//
// Phase-1 stub CLI for the `tise_solver` executable described in
// docs/SDD.md §7.2.1 (Controller to TISE Solver). This file implements
// ONLY the subprocess contract -- argument parsing, config loading, output
// directory creation, and the expected output file set -- with placeholder
// file contents. Real TISE physics/numerics are deferred to §10.2 Phase 4
// and will replace this file's internals; nothing here should be treated
// as a foundation to build real numerics on top of.
//
// Usage:
//   tise_solver --config <config.yaml> --output-dir <output_dir>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <yaml-cpp/yaml.h>

namespace
{

constexpr const char *USAGE = "usage: tise_solver --config <path> --output-dir <path>";

// Placeholder dimensions for the matrix-shaped stub outputs (see
// writeGridFile() below), kept as named constants so the made-up dimensions
// are consistent across files instead of each being independently guessed.
// These are shape-only choices for a Phase-1 contract stub -- they do NOT
// encode any real LAPACK banded-storage convention (undocumented until a
// later phase; see docs/SDD.md §10.2 Phase 4) and are not physically
// meaningful basis sizes.
constexpr int kPlaceholderNumEigenstates = 5; // eigenvalues.dat rows; eigenvectors.dat cols
constexpr int kPlaceholderBasisSize = 8;      // eigenvectors.dat rows; hamiltonian.dat/overlap.dat cols
constexpr int kPlaceholderBandwidth = 3;      // hamiltonian.dat/overlap.dat rows (< basis size: suggests banded/compact storage)

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

// Write a 2-column (index, value) placeholder file: `numRows` rows of
// "<index>  0.0" beneath a header comment line. This is the original
// writeStubDataFile() body, parameterized by row count instead of hardcoded
// to 2. Used for outputs that docs/SDD.md §6.3 documents as genuinely
// index/value pairs: eigenvalues.dat (index, E_n) and
// continuum_state_NNN.dat (x, psi_eps(x)). Content is not physically
// meaningful -- this is a Phase-1 contract stub, not real physics (see
// docs/SDD.md §10.2 Phase 4).
void writeIndexValueFile(const std::filesystem::path &path, const std::string &description, int numRows)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("failed to open '" + path.string() + "' for writing");

    out << "# " << description << " -- Phase 1 stub, no physics yet\n";
    for (int i = 0; i < numRows; ++i)
        out << i << "  0.0\n";
}

// Write a `numRows` x `numCols` rectangular grid of placeholder 0.0 values,
// whitespace-separated with no index column, beneath the same header-comment
// style as writeIndexValueFile(). Used for the files §6.3 documents as
// matrix/banded-matrix content -- eigenvectors.dat, hamiltonian.dat,
// overlap.dat -- so their on-disk shape is structurally a matrix rather than
// an index-value list. This deliberately does NOT encode any real LAPACK
// banded-storage convention (bandwidth formula, compact-storage layout) --
// that's undocumented until a later phase (real TISE numerics); only the
// row/column shape is meaningful here, the placeholder values stay
// meaningless zeros.
void writeGridFile(const std::filesystem::path &path, const std::string &description, int numRows, int numCols)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("failed to open '" + path.string() + "' for writing");

    out << "# " << description << " -- Phase 1 stub, no physics yet\n";
    for (int r = 0; r < numRows; ++r)
    {
        for (int c = 0; c < numCols; ++c)
            out << (c == 0 ? "" : "  ") << "0.0";
        out << "\n";
    }
}

// Write phase_shifts.dat's 3-column (eps_i, delta(eps_i), d(delta)/dE)
// placeholder: `numRows` rows of "0.0  0.0  0.0" beneath the same
// header-comment style. `numRows` must be the same n_energies value that
// drives the continuum_state_NNN.dat file count in writeContinuumOutputs()
// below -- §6.3 documents one phase_shifts.dat row per energy point,
// mirroring one continuum_state file per energy point.
void writePhaseShiftsFile(const std::filesystem::path &path, int numRows)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("failed to open '" + path.string() + "' for writing");

    out << "# phase_shifts.dat: eps_i, delta(eps_i), d(delta)/dE -- Phase 1 stub, no physics yet\n";
    for (int i = 0; i < numRows; ++i)
        out << "0.0  0.0  0.0\n";
}

// Zero-pad `value` to at least `width` digits, e.g. zeroPadded(1, 3) == "001".
std::string zeroPadded(int value, int width)
{
    std::ostringstream oss;
    oss << std::setw(width) << std::setfill('0') << value;
    return oss.str();
}

// Write the 5 files the Controller<->TISE contract (docs/SDD.md §6.3,
// §7.2.1) always expects, regardless of continuum settings.
void writeCoreOutputs(const std::filesystem::path &outputDir)
{
    writeIndexValueFile(outputDir / "eigenvalues.dat", "eigenvalues.dat: index, E_n", kPlaceholderNumEigenstates);
    writeGridFile(outputDir / "eigenvectors.dat", "eigenvectors.dat: columns are c_n coefficient vectors",
                  kPlaceholderBasisSize, kPlaceholderNumEigenstates);
    writeGridFile(outputDir / "hamiltonian.dat", "hamiltonian.dat: H matrix (banded)",
                  kPlaceholderBandwidth, kPlaceholderBasisSize);
    writeGridFile(outputDir / "overlap.dat", "overlap.dat: S matrix (banded)",
                  kPlaceholderBandwidth, kPlaceholderBasisSize);

    const std::filesystem::path warningsPath = outputDir / "warnings.json";
    std::ofstream warnings(warningsPath);
    if (!warnings)
        throw std::runtime_error("failed to open '" + warningsPath.string() + "' for writing");
    // No real physics warnings exist until a later phase; an empty JSON
    // array is the well-formed placeholder a Python Controller can safely
    // json.load() today, ready to be populated once warnings exist.
    warnings << "[]\n";
}

// Validate and extract config["tise"]["continuum"]["n_energies"]. Deliberately
// separated from writeContinuumOutputs() and called BEFORE writeCoreOutputs()
// in main(): the Controller<->TISE contract (docs/SDD.md §7.2.1) treats any
// non-zero exit as a "hard failure" that must leave no partial output behind,
// so everything that could invalidate this run has to be validated before the
// first byte of any output file is written.
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

// Write phase_shifts.dat plus one continuum_state_NNN.dat per energy point.
// Only called when config["tise"]["continuum"]["enabled"] is true, with an
// nEnergies already validated by validateNEnergies().
void writeContinuumOutputs(const std::filesystem::path &outputDir, int nEnergies)
{
    writePhaseShiftsFile(outputDir / "phase_shifts.dat", nEnergies);

    for (int i = 1; i <= nEnergies; ++i)
    {
        const std::string filename = "continuum_state_" + zeroPadded(i, 3) + ".dat";
        writeIndexValueFile(outputDir / filename, filename + ": x, psi_eps(x)", 2);
    }
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

        // yaml-cpp idiom: an undefined/missing Node converts to false
        // without throwing, so configs that omit tise.continuum entirely,
        // or set enabled: false, correctly skip this branch.
        const bool continuumEnabled = config["tise"] && config["tise"]["continuum"] &&
                                       config["tise"]["continuum"]["enabled"] &&
                                       config["tise"]["continuum"]["enabled"].as<bool>();

        // Validate everything this run will need BEFORE writing anything:
        // if continuum is enabled but n_energies is missing/invalid, fail
        // here, before writeCoreOutputs() has written a single file (see
        // docs/SDD.md §7.2.1: a non-zero exit is a "hard failure" and must
        // not leave partial output on disk).
        int nEnergies = 0;
        if (continuumEnabled)
            nEnergies = validateNEnergies(config["tise"]["continuum"]);

        writeCoreOutputs(outputDir);

        if (continuumEnabled)
            writeContinuumOutputs(outputDir, nEnergies);
    }
    catch (const std::exception &e)
    {
        std::cerr << "tise_solver: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
