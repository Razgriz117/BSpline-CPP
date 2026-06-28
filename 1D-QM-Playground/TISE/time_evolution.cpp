#include "time_evolution.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace tevol
{

namespace
{
    constexpr double PI = 3.14159265358979323846;
}

double gaussianWavepacket(double x, double r0, double mOmega, double hbar)
{
    const double dx = x - r0;
    return std::pow(mOmega / (PI * hbar), 0.25) *
           std::exp(-mOmega * dx * dx / (2.0 * hbar));
}

Eigen::VectorXd computeGaussianOverlaps(const bspline::BSpline &bs,
                                         int nEn,
                                         int order,
                                         double r0,
                                         double mOmega,
                                         double hbar)
{
    Eigen::VectorXd B_G(nEn);

    bspline::D2DFun fGauss = [r0, mOmega, hbar](double x, const double *) {
        return gaussianWavepacket(x, r0, mOmega, hbar);
    };

    for (int iBs2 = 2; iBs2 <= nEn + 1; ++iBs2)
    {
        double sum = 0.0;
        int iBs1Min = std::max(2, iBs2 - order + 1);
        for (int iBs1 = iBs1Min; iBs1 <= iBs2; ++iBs1)
            sum += bs.integral(fGauss, iBs2, iBs1);

        B_G(iBs2 - 2) = sum; // 0-based: iBs2=2 -> index 0
    }

    return B_G;
}

Eigen::VectorXd projectToEigenBasis(const Eigen::MatrixXd &Cinv,
                                     const Eigen::VectorXd &B_G)
{
    return Cinv * B_G;
}

Eigen::VectorXcd timeEvolveState(const Eigen::VectorXd &Phi_G,
                                  const std::vector<Real> &eval,
                                  double t,
                                  double hbar)
{
    const int n = static_cast<int>(eval.size());
    Eigen::VectorXcd V_t(n);
    for (int k = 0; k < n; ++k)
    {
        std::complex<double> phase(0.0, -eval[k] * t / hbar);
        V_t(k) = std::exp(phase) * Phi_G(k);
    }
    return V_t;
}

Eigen::VectorXcd transformToSpaceBasis(const Eigen::MatrixXd &C,
                                        const Eigen::VectorXcd &V_t)
{
    return C * V_t;
}

void writeTimestep(std::ostream &out,
                   const bspline::BSpline &bs,
                   const Eigen::VectorXcd &CV_padded,
                   int nBSplines,
                   int npts,
                   Real rMin,
                   Real rMax)
{
    out << std::scientific << std::setprecision(16);

    std::vector<Real> real_parts(nBSplines), imag_parts(nBSplines);
    for (int i = 0; i < nBSplines; ++i)
    {
        real_parts[i] = CV_padded(i).real();
        imag_parts[i] = CV_padded(i).imag();
    }

    for (int ix = 1; ix <= npts; ++ix)
    {
        Real x = rMin + (rMax - rMin) *
                 static_cast<Real>(ix - 1) / static_cast<Real>(npts - 1);
        double re = bs.eval(x, real_parts.data(), nBSplines);
        double im = bs.eval(x, imag_parts.data(), nBSplines);
        out << " " << std::setw(24) << x
            << " " << std::setw(24) << re
            << " " << std::setw(24) << im << "\n";
    }
}

void runTimeEvolution(const bspline::BSpline &bs,
                      const tise::EigenResult &er,
                      int nBSplines,
                      Real rMin,
                      Real rMax,
                      int timeSteps,
                      double dt,
                      double r0,
                      double mOmega,
                      double hbar,
                      const std::string &outputDir)
{
    const int nEn = er.dim;
    const int order = bs.getOrder();

    Eigen::VectorXd B_G = computeGaussianOverlaps(bs, nEn, order, r0, mOmega, hbar);

    Eigen::Map<const Eigen::MatrixXd> C(er.vectors.data(), nEn, nEn);
    Eigen::MatrixXd Cinv = C.inverse();

    Eigen::VectorXd Phi_G = projectToEigenBasis(Cinv, B_G);

    for (int step = 0; step < timeSteps; ++step)
    {
        double t = step * dt;

        Eigen::VectorXcd V_t  = timeEvolveState(Phi_G, er.values, t, hbar);
        Eigen::VectorXcd CV_t = transformToSpaceBasis(C, V_t);

        Eigen::VectorXcd CV_padded(nBSplines);
        CV_padded.setZero();
        CV_padded.segment(1, nEn) = CV_t;

        std::ostringstream oss;
        oss << outputDir << "/Timestep_"
            << std::setw(3) << std::setfill('0') << step;

        std::ofstream out(oss.str());
        if (!out)
            throw std::runtime_error("Cannot open " + oss.str());

        writeTimestep(out, bs, CV_padded, nBSplines, 301, rMin, rMax);
    }
}

} // namespace tevol
