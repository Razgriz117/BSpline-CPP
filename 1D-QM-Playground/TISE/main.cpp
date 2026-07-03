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
constexpr int L = 1;

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

std::map<std::string, std::string> parsePiecewise(const std::string& arg) {
    nlohmann::json j = nlohmann::json::parse(arg);          // arg must be valid JSON array
    std::map<std::string, std::string> domainToFunction;

    for (const auto& piece : j) {
        domainToFunction[piece.at("domain").get<std::string>()] =
            piece.at("function").get<std::string>();
    }
    return domainToFunction;
}

int main(int argc, char *argv[])
{
    // Parse the potential input string into a map
    auto potential = parsePiecewise(argv[1]);

    std::cout << "Potential is: " << std::endl;
    for (const auto& [domain, fn] : potential) {
        std::cout << "\t" << domain << ": " << fn << "\n";
    }

    // ------------------------------------------------------------------
    // Project Part 1: solve the Time-Independent Schrödinger Equation
    // ------------------------------------------------------------------
    tise::EigenResult er;
    try
    {
        er = tise::solveTISE(BS_NNODS, BS_ORDER, BS_GRMIN, BS_GRMAX, L, potential);
    }
    catch (const std::exception &e)
    {
        std::cerr << "TISE solver failed: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // Rebuild the B-spline basis (needed for eigenfunction evaluation and Part 2)
    auto grid = tise::buildUniformRadialGrid(BS_NNODS, BS_GRMIN, BS_GRMAX);
    bspline::BSpline bs;
    if (bs.init(BS_NNODS, BS_ORDER, grid) != 0)
    {
        std::cerr << "BSpline::init failed\n";
        return EXIT_FAILURE;
    }
    const int nBSplines = bs.getNBSplines();
    const int nEn       = er.dim;

    // Print accurate eigenvalues and write eigenstate files
    std::cout << std::scientific << std::setprecision(16);
    int iEn;
    for (iEn = 1; iEn <= nEn; ++iEn)
    {
        double eig = er.values[iEn - 1];
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
        auto coeffs = tise::eigenstateCoefficients(er.vectors, iEn, nEn, nBSplines);
        tise::writeEigenstate(out, bs, coeffs, NPTS_EIGENSTATE, BS_GRMIN, BS_GRMAX);
    }
    std::cout << "Number of Accurate Eigenvalues : " << (iEn - 1) << "\n";

    // ------------------------------------------------------------------
    // Project Part 2: propagate a Gaussian wavepacket in time
    // ------------------------------------------------------------------
    try
    {
        tevol::runTimeEvolution(bs, er, nBSplines,
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
