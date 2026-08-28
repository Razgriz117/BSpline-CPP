#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include <nlohmann/json.hpp>
#include "BSpline.hpp"
#include "tise.hpp"
#include "time_evolution.hpp"

#ifndef TIMESTEPS_DIR
     #define TIMESTEPS_DIR "./timesteps"
#endif

// Angular momentum quantum number (change to 0, 1, 2, ...)
constexpr int L = 0;

// Physical parameters for the Gaussian wavepacket (atomic units)
constexpr double HBAR             = 1.0;
constexpr double PARTICLE_MASS    = 3.75;
constexpr double INITIAL_POSITION = 10.0;
constexpr double K_SPRING         = 1.0 / 8.0;
const     double M_OMEGA          = std::sqrt(PARTICLE_MASS * K_SPRING);

// B-spline basis parameters
constexpr int    BS_NNODS = 51;
constexpr int    BS_ORDER = 12;
constexpr double BS_GRMIN = 0.0;
constexpr double BS_GRMAX = 100.0;

// Output parameters
constexpr double ERROR_THRESHOLD = 1.0e-10;
constexpr int    NPTS_EIGENSTATE = 301;
constexpr int    TIME_STEPS      = 1000;
constexpr double DT              = 0.3;

// Continuum energy-grid parameters (REQ-F-040); see buildEnergyGrid.
constexpr double E_THRESHOLD = 0.0;
constexpr double E_MAX       = 2.0;
constexpr int    N_E         = 50;

// Parse a JSON array of {"domain": ..., "function": ...} objects (as produced by
// the controller module from a YAML "potential" list) into a domain-string -> expression-string
// map. Each entry describes one piece of a piecewise-defined potential; the
// expression is evaluated later by muparser in tise::evaluateFunction.
std::map<std::string, std::string> parsePiecewise(const std::string& arg) {
    nlohmann::json function_array = nlohmann::json::parse(arg);  // arg: valid JSON array
    std::map<std::string, std::string> domainToFunction;

    for (const auto &piece : function_array)
    {
        domainToFunction[piece.at("domain").get<std::string>()] =
            piece.at("function").get<std::string>();
    }
    return domainToFunction;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " '<JSON array of {\"domain\":...,\"function\":...} pieces>'\n";
        return EXIT_FAILURE;
    }

    // argv[1]: JSON-encoded array describing the piecewise potential, e.g.
    // [{"domain": "[0, 20)", "function": "x"}, {"domain": "[20, 40]", "function": "x^2"}]
    auto potential = parsePiecewise(argv[1]);

    std::cout << "Potential is: " << std::endl;
    for (const auto& [domain, fn] : potential) {
        std::cout << "\t" << domain << ": " << fn << "\n";
    }

    // ------------------------------------------------------------------
    // Project Part 1: solve the Time-Independent Schrödinger Equation
    // ------------------------------------------------------------------
    tise::SolveTISEResult sol;
    try
    {
        sol = tise::solveTISE(BS_NNODS, BS_ORDER, BS_GRMIN, BS_GRMAX, L, potential,
                              E_THRESHOLD, E_MAX, N_E);
    }
    catch (const std::exception &e)
    {
        std::cerr << "TISE solver failed: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // solveTISE already built (and diagonalized against) the basis/drop-set
    // below -- reuse directly instead of independently rebuilding, which
    // would silently diverge now that solveTISE may build a strategic grid
    // and/or a non-classic drop-set (engineer-a-plan-A4-wiring.md Gaps 3/4).
    const bspline::BSpline &bs = sol.bs;
    const int nBSplines = sol.nBSplines;
    const int nEn       = sol.eigen.dim;

    // Print accurate eigenvalues and write eigenstate files
    std::cout << std::scientific << std::setprecision(16);
    int iEn;
    for (iEn = 1; iEn <= nEn; ++iEn)
    {
        double eig = sol.eigen.values[iEn - 1];
        double err = tise::eigenvalueError(eig, iEn, L);
        if (err > ERROR_THRESHOLD) break;

        std::cout << std::setw(4) << iEn << "  "
                  << std::setw(24) << eig << "  "
                  << std::setw(24) << err << "\n";

        std::ostringstream oss;
        oss << "EigenState_" << std::setw(3) << std::setfill('0') << iEn;
        std::ofstream out(oss.str());
        if (!out)
        {
            std::cerr << "Cannot open " << oss.str() << "\n";
            return EXIT_FAILURE;
        }
        auto coeffs = tise::eigenstateCoefficients(sol.eigen.vectors, iEn, nEn, nBSplines, sol.dropSet);
        tise::writeEigenstate(out, bs, coeffs, NPTS_EIGENSTATE, BS_GRMIN, BS_GRMAX);
    }
    std::cout << "Number of Accurate Eigenvalues : " << (iEn - 1) << "\n";

    // ------------------------------------------------------------------
    // Project Part 2: propagate a Gaussian wavepacket in time
    // ------------------------------------------------------------------
    // NOTE (known limitation -- tise-task-breakdown.md Sec. 4 item 4): this call
    // is unconditional -- time evolution always runs after the TISE solve above,
    // regardless of intent. config.yaml already defines a matching run.run_tdse
    // flag (docs/SDD.md Sec. 6.1), but main.cpp has no YAML parsing of its own and
    // never reads it; gating this call on that flag needs the (out-of-scope)
    // config-driven Controller<->TISE plumbing from the interface phases. Flagged
    // here, not fixed, per the cleanup task's own scope.
    try
    {
        tevol::runTimeEvolution(sol.bs, sol.eigen, sol.nBSplines,
                                BS_GRMIN, BS_GRMAX,
                                TIME_STEPS, DT,
                                INITIAL_POSITION, M_OMEGA, HBAR,
                                TIMESTEPS_DIR);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Time evolution failed: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return 0;
}
