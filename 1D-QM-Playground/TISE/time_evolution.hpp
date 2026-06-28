#pragma once

#include <complex>
#include <iosfwd>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "BSpline.hpp"
#include "tise.hpp"

namespace tevol
{

using Real = bspline::Real;

// Ground state of a harmonic oscillator (normalised Gaussian wavepacket).
// mOmega = m * omega (in a.u.); hbar = 1.0 in atomic units.
double gaussianWavepacket(double x, double r0, double mOmega, double hbar);

// Compute <B_i | G> for i = 2..nEn+1 (the active B-spline range).
// Returns a dense vector of length nEn — no padding.
Eigen::VectorXd computeGaussianOverlaps(const bspline::BSpline &bs,
                                         int nEn,
                                         int order,
                                         double r0,
                                         double mOmega,
                                         double hbar);

// Project B-spline overlaps to the eigenstate basis: Phi_G = Cinv * B_G.
Eigen::VectorXd projectToEigenBasis(const Eigen::MatrixXd &Cinv,
                                     const Eigen::VectorXd &B_G);

// Apply diagonal time-evolution: V(t)[k] = exp(-i E_k t / hbar) * Phi_G[k].
Eigen::VectorXcd timeEvolveState(const Eigen::VectorXd &Phi_G,
                                  const std::vector<Real> &eval,
                                  double t,
                                  double hbar);

// Transform from eigenstate basis back to B-spline basis: CV_t = C * V_t.
Eigen::VectorXcd transformToSpaceBasis(const Eigen::MatrixXd &C,
                                        const Eigen::VectorXcd &V_t);

// Write one timestep to `out`: npts lines of "x  Re(psi)  Im(psi)".
// CV_padded has length nBSplines (including boundary zeros).
void writeTimestep(std::ostream &out,
                   const bspline::BSpline &bs,
                   const Eigen::VectorXcd &CV_padded,
                   int nBSplines,
                   int npts,
                   Real rMin,
                   Real rMax);

// Top-level time evolution: propagate a Gaussian wavepacket for timeSteps
// steps of size dt and write each frame to outputDir/Timestep_NNN.
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
                      const std::string &outputDir);

} // namespace tevol
