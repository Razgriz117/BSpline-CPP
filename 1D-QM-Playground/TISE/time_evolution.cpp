#include "time_evolution.hpp"

#include <algorithm>
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

    // <B_i|G> is obtained from the partition-of-unity identity sum_j B_j(x) = 1:
    //
    //   <B_i|G> = int B_i(x) G(x) dx
    //           = int B_i(x) G(x) [sum_j B_j(x)] dx
    //           = sum_j integral(G, i, j)
    //
    // BSpline::integral only ever forms TWO-spline matrix elements, so this
    // identity is the route to a single-spline overlap. The sum must run over
    // EVERY j whose support meets B_i's -- the full band |i - j| <= order-1,
    // clamped to the physical range [1, nBSplines] -- otherwise the partition
    // of unity is truncated and the identity does not hold.
    //
    // This loop previously reused the BANDED-MATRIX FILL bounds from
    // Template.f90/main.cpp (j = max(2, i-order+1) .. i), which walk only the
    // lower half of the band because a symmetric matrix needs nothing more.
    // That is correct for filling H and S; it is wrong here, and undercounted
    // <B_i|G> by a factor of ~3 precisely where the wavepacket carries its
    // amplitude. Guarded by ComputeGaussianOverlapsTest.MatchesDirectQuadrature.
    const int nBSplines = bs.getNBSplines();

    for (int iBs2 = 2; iBs2 <= nEn + 1; ++iBs2)
    {
        double sum = 0.0;
        const int iBs1Min = std::max(1, iBs2 - order + 1);
        const int iBs1Max = std::min(nBSplines, iBs2 + order - 1);
        for (int iBs1 = iBs1Min; iBs1 <= iBs1Max; ++iBs1)
            sum += bs.integral(fGauss, iBs2, iBs1);

        B_G(iBs2 - 2) = sum; // 0-based: iBs2=2 -> index 0
    }

    return B_G;
}

Eigen::VectorXd projectToEigenBasis(const Eigen::MatrixXd &M,
                                     const Eigen::VectorXd &B_G)
{
    return M * B_G;
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

    // Spectral amplitudes a_n = <phi_n|G>, for |phi_n> = sum_i C(i,n) |B_i>.
    //
    // DSBGV returns S-ORTHONORMAL eigenvectors (C^T S C = I), so the phi_n are
    // already orthonormal under the B-spline overlap matrix and the amplitude
    // is just the C-weighted sum of the B-spline overlaps in B_G:
    //
    //   a_n = <phi_n|G> = sum_i C(i,n) <B_i|G>   =>   Phi_G = C^T B_G
    //
    // NOT C^{-1} B_G. S-orthonormality means C^{-1} = C^T S, so routing through
    // the inverse applies one spurious extra factor of S. (C^{-1} is the right
    // operator for a B-spline COEFFICIENT vector c -- and B_G could be turned
    // into one via c = S^{-1} B_G -- but B_G holds overlaps, not coefficients,
    // and the two S factors then cancel to leave exactly C^T.) Transposing is
    // also O(n^2) against C.inverse()'s O(n^3), and avoids inverting a matrix
    // that never needed inverting.
    Eigen::VectorXd Phi_G = projectToEigenBasis(C.transpose(), B_G);

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
