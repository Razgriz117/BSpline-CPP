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

// Write a minimal placeholder data file: a header comment line plus a couple
// of placeholder rows. Content is not physically meaningful -- this is a
// Phase-1 contract stub, not real physics (see docs/SDD.md §10.2 Phase 4).
void writeStubDataFile(const std::filesystem::path &path, const std::string &description)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("failed to open '" + path.string() + "' for writing");

    out << "# " << description << " -- Phase 1 stub, no physics yet\n";
    out << "0  0.0\n";
    out << "1  0.0\n";
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
    writeStubDataFile(outputDir / "eigenvalues.dat", "eigenvalues.dat: index, E_n");
    writeStubDataFile(outputDir / "eigenvectors.dat", "eigenvectors.dat: columns are c_n coefficient vectors");
    writeStubDataFile(outputDir / "hamiltonian.dat", "hamiltonian.dat: H matrix (banded)");
    writeStubDataFile(outputDir / "overlap.dat", "overlap.dat: S matrix (banded)");

    const std::filesystem::path warningsPath = outputDir / "warnings.json";
    std::ofstream warnings(warningsPath);
    if (!warnings)
        throw std::runtime_error("failed to open '" + warningsPath.string() + "' for writing");
    // No real physics warnings exist until a later phase; an empty JSON
    // array is the well-formed placeholder a Python Controller can safely
    // json.load() today, ready to be populated once warnings exist.
    warnings << "[]\n";
}

// Write phase_shifts.dat plus one continuum_state_NNN.dat per energy point.
// Only called when config["tise"]["continuum"]["enabled"] is true.
void writeContinuumOutputs(const std::filesystem::path &outputDir, const YAML::Node &continuumNode)
{
    // Missing or wrong-typed n_energies throws YAML::TypedBadConversion<int>
    // (or YAML::InvalidNode if the key is absent entirely) -- both derive
    // from std::exception and propagate to main()'s catch block unmodified.
    const int nEnergies = continuumNode["n_energies"].as<int>();
    if (nEnergies <= 0)
        throw std::runtime_error("tise.continuum.n_energies must be a positive integer");

    writeStubDataFile(outputDir / "phase_shifts.dat", "phase_shifts.dat: eps_i, delta(eps_i), d(delta)/dE");

    for (int i = 1; i <= nEnergies; ++i)
    {
        const std::string filename = "continuum_state_" + zeroPadded(i, 3) + ".dat";
        writeStubDataFile(outputDir / filename, filename + ": x, psi_eps(x)");
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

        writeCoreOutputs(outputDir);

        // yaml-cpp idiom: an undefined/missing Node converts to false
        // without throwing, so configs that omit tise.continuum entirely,
        // or set enabled: false, correctly skip this branch.
        const bool continuumEnabled = config["tise"] && config["tise"]["continuum"] &&
                                       config["tise"]["continuum"]["enabled"] &&
                                       config["tise"]["continuum"]["enabled"].as<bool>();

        if (continuumEnabled)
            writeContinuumOutputs(outputDir, config["tise"]["continuum"]);
    }
    catch (const std::exception &e)
    {
        std::cerr << "tise_solver: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
