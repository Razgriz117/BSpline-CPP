// This file is a C++ translation and adaptation of the original Fortran
// ModuleBSpline code by Luca Argenti.
//
// ---------------------------------------------------------------------
// Original copyright and license notice:
//
// ModuleBSplines, Copyright (C) 2020 Luca Argenti, PhD - Some Rights Reserved
// ModuleBSplines is licensed under a
// Creative Commons Attribution-ShareAlike 4.0 International License.
// A copy of the license is available at <http://creativecommons.org/licenses/by-nd/4.0/>.
//
// Luca Argenti is Associate Professor of Physics, Optics and Photonics
// Department of Physics and the College of Optics
// University of Central Florida
// 4111 Libra Drive, Orlando, Florida 32816, USA
// email: luca.argenti@ucf.edu
//
// ---------------------------------------------------------------------
//
// C++ translation and modifications:
//   Copyright (C) 2025  Martin de Salterain
//   Copyright (C) 2025  Hui Xian Grace Lim
//
// This C++ version is a derivative work based on ModuleBSplines and is
// distributed under the same Creative Commons license as the original.
// Please refer to the original license terms for details.

#include "BSpline.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <vector>

namespace bspline
{

     // ========================
     // Module-level constants
     // ========================

     namespace
     {

          constexpr int NGAUSS = 64;
          constexpr Real PI = 3.14159265358979323846;
          constexpr Real GAUSS_POINTS_CONV_THRESH = 2.0e-15;
          constexpr int DIMFACT = 127;
          constexpr Real NOD_THRESHOLD = 1.0e-15;

          // Gaussian quadrature nodes and weights on [0,1]
          Real gaussPoints[NGAUSS];
          Real gaussWeights[NGAUSS];

          // Factorials and binomials (for completeness, as in Fortran).
          Real fact[DIMFACT + 1];
          constexpr int TAR_SIZE = ((DIMFACT + 1) * (DIMFACT + 2)) / 2;
          Real tar[TAR_SIZE];

          // Partial factorial matrix parfactMat(0:DIMFACT, 0:DIMFACT)
          Real parfactMat[(DIMFACT + 1) * (DIMFACT + 1)];

          inline Real &parfact(int row, int col)
          {
               // row, col in 0..DIMFACT (row-major storage)
               return parfactMat[row * (DIMFACT + 1) + col];
          }

          inline Real partialFactorial(int k, int n)
          {
               // Fortran uses parfactmat(k,n)
               // We stored parfactMat(row=n, col=k) as parfact(n,k)
               // so here we must return parfact(k,n).
               return parfact(k, n);
          }

          bool moduleInitialized = false;

          // Initialize Gauss-Legendre nodes/weights on [0,1], mirroring InitGaussPoints().
          void initGaussPoints()
          {
               for (int i = 1; i <= (NGAUSS + 1) / 2; ++i)
               {
                    // Initial guess for root on [-1,1]
                    Real z = std::cos(PI * (static_cast<Real>(i) - 0.25) / (static_cast<Real>(NGAUSS) + 0.5));
                    Real z1, p1, p2, p3, pp;

                    // Newton-Raphson
                    while (true)
                    {
                         p1 = 1.0;
                         p2 = 0.0;
                         for (int j = 1; j <= NGAUSS; ++j)
                         {
                              p3 = p2;
                              p2 = p1;
                              p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
                         }
                         pp = NGAUSS * (z * p1 - p2) / (z * z - 1.0);
                         z1 = z;
                         z = z1 - p1 / pp;
                         if (std::abs(z - z1) <= GAUSS_POINTS_CONV_THRESH)
                         {
                              break;
                         }
                    }

                    // Map [-1,1] -> [0,1]
                    Real x1 = (1.0 - z) * 0.5;
                    Real x2 = (1.0 + z) * 0.5;
                    Real w = 1.0 / ((1.0 - z * z) * pp * pp);

                    gaussPoints[i - 1] = x1;
                    gaussPoints[NGAUSS - i] = x2;
                    gaussWeights[i - 1] = w;
                    gaussWeights[NGAUSS - i] = w;
               }
          }

          // Initialize factorials, binomials, parfactMat(), mirroring InitFactorials().
          void initFactorials()
          {
               // Factorials
               fact[0] = 1.0;
               for (int i = 1; i <= DIMFACT; ++i)
               {
                    fact[i] = fact[i - 1] * static_cast<Real>(i);
               }

               // Binomial coefficients tar[]
               int k = 0;
               for (int i = 0; i <= DIMFACT; ++i)
               {
                    for (int j = 0; j <= i; ++j)
                    {
                         tar[k] = fact[i] / (fact[i - j] * fact[j]);
                         ++k;
                    }
               }

               // Initialize parfactMat to 1.0
               for (int i = 0; i < (DIMFACT + 1) * (DIMFACT + 1); ++i)
               {
                    parfactMat[i] = 1.0;
               }

               // parfactMat(n,k) = Π_{i=n-k+1..n} i (n=1..DIMFACT, k=1..n)
               for (int n = 1; n <= DIMFACT; ++n)
               {
                    for (int kk = 1; kk <= n; ++kk)
                    {
                         Real val = 1.0;
                         for (int i = n - kk + 1; i <= n; ++i)
                         {
                              val *= static_cast<Real>(i);
                         }
                         parfact(n, kk) = val;
                    }
               }
          }

          void initSplineModule()
          {
               if (!moduleInitialized)
               {
                    initGaussPoints();
                    initFactorials();
                    moduleInitialized = true;
               }
          }

          inline Real upperLimitTo(Real x)
          {
               const Real eps = std::numeric_limits<Real>::epsilon();
               const Real tin = std::numeric_limits<Real>::min() * (1.0 + eps);
               return (x + tin) * (1.0 + eps);
          }

          inline Real lowerLimitTo(Real x)
          {
               const Real eps = std::numeric_limits<Real>::epsilon();
               const Real tin = std::numeric_limits<Real>::min() * (1.0 + eps);
               return (x - tin) * (1.0 - eps);
          }

          // Static saved index used in whichInterval() (mirrors Fortran "save :: i1").
          int saved_i1 = 0;

     } // anonymous namespace

     // Internal unity function used for normalization within BSpline::init().
     Real unityFunction(Real /*x*/, const Real * /*parvec*/)
     {
          return 1.0;
     }

     // ========================
     // BSpline member functions
     // ========================

     int BSpline::init(int numberOfNodes, int order, const std::vector<Real> &gridIn)
     {
          initSplineModule();

          // Check basic parameters (like CheckNodesAndOrder in Fortran).
          if (numberOfNodes < 2)
               return -1;
          if (order < 1)
               return -2;

          nNodes_ = numberOfNodes;
          order_ = order;
          nBSplines_ = numberOfNodes + order - 2;

          if (static_cast<int>(gridIn.size()) < numberOfNodes)
          {
               return -4; // grid size insufficient
          }

          // Ensure grid is non-decreasing.
          for (int i = 1; i < numberOfNodes; ++i)
          {
               if (gridIn[i] < gridIn[i - 1])
               {
                    return 1; // Fortran's error when grid(i) < grid(i-1)
               }
          }

          // Extended grid: indices from -(order_-1) .. (nNodes_-1 + order_-1).
          grid_.assign(nNodes_ + 2 * order_ - 2, 0.0);

          // grid(1-order : -1) = gridIn(1)
          for (int i = 1 - order_; i <= -1; ++i)
          {
               gridAt(i) = gridIn[0];
          }
          // grid(0 : nNodes_-1) = distinct nodes
          for (int i = 0; i <= nNodes_ - 1; ++i)
          {
               gridAt(i) = gridIn[i];
          }
          // grid(nNodes_ : nNodes_ + order_ - 2) = gridIn(nNodes_)
          for (int i = nNodes_; i <= nNodes_ + order_ - 2; ++i)
          {
               gridAt(i) = gridIn[nNodes_ - 1];
          }

          // Allocate coefficients c(0:order_-1, 0:order_-1, 1:nBSplines_).
          coeffs_.assign(order_ * order_ * nBSplines_, 0.0);

          std::vector<Real> singleCoeffs(order_ * order_);

          for (int Bs = 1; Bs <= nBSplines_; ++Bs)
          {
               const Real *vec = &grid_[idxGrid(Bs - order_)]; // vec(1) in Fortran
               computeCoefficientsSingleBSpline(order_, vec, singleCoeffs.data());
               for (int m = 0; m < order_; ++m)
               {
                    for (int j = 0; j < order_; ++j)
                    {
                         coeff(m, j, Bs) = singleCoeffs[m + order_ * j];
                    }
               }
          }

          // At this point, B-spline shape is fully defined, so we can safely
          // mark the object as initialized **before** calling integral() for
          // normalization (this was the source of your crash).
          initialized_ = true;

          // Compute L^2 normalization factors f(i) = 1/sqrt(∫ B_i^2 dx).
          normalization_.assign(nBSplines_, 0.0);
          D2DFun funUnity = unityFunction;
          for (int i = 1; i <= nBSplines_; ++i)
          {
               Real I = integral(funUnity, i, i);
               normalization_[i - 1] = 1.0 / std::sqrt(I);
          }

          return 0;
     }

     void BSpline::free()
     {
          grid_.clear();
          normalization_.clear();
          coeffs_.clear();
          nNodes_ = 0;
          order_ = 0;
          nBSplines_ = 0;
          initialized_ = false;
     }

     void BSpline::checkInitialization() const
     {
          if (!initialized_)
          {
               throw std::runtime_error("BSpline used before initialization");
          }
     }

     Real BSpline::getNodePosition(int node) const
     {
          checkInitialization();

          if (node <= 1)
          {
               // Fortran: Position = SplineSet%Grid(0)
               return gridAt(0);
          }
          else if (node >= nNodes_)
          {
               // Fortran: Position = SplineSet%Grid(NNodes-1)
               return gridAt(nNodes_ - 1);
          }
          else
          {
               // Fortran: Position = SplineSet%Grid(Node-1)
               return gridAt(node - 1);
          }
     }

     Real BSpline::getNormFactor(int bsIndex) const
     {
          checkInitialization();

          if (bsIndex < 1 || bsIndex > nBSplines_)
          {
               // Fortran BSplineNormalizationFactor returns 0 if out-of-range
               return Real(0);
          }
          // Fortran indexing: f(Bs) -> normalization_[Bs-1]
          return normalization_[bsIndex - 1];
     }

     void BSpline::computeCoefficientsSingleBSpline(int order,
                                                    const Real *vec,
                                                    Real *coefficients) const
     {
          const int dim = order;

          // coe(m, inter, k, tempBs), m=0..order-1, inter=0..order-1,
          // k=1..order, tempBs=1..order
          std::vector<Real> coe(dim * dim * dim * dim, 0.0);

          auto idxCoe = [dim](int m, int inter, int k, int tempBs)
          {
               // 4D index flattening (Fortran-like):
               // m + dim*(inter + dim*((k-1) + dim*(tempBs-1)))
               return m + dim * (inter + dim * ((k - 1) + dim * (tempBs - 1)));
          };

          // Initial condition: coe(0,0,1,:) = 1
          for (int tempBs = 1; tempBs <= order; ++tempBs)
          {
               coe[idxCoe(0, 0, 1, tempBs)] = 1.0;
          }

          // k = 2..order
          for (int kOrder = 2; kOrder <= order; ++kOrder)
          {
               for (int tempBs = 1; tempBs <= order - kOrder + 1; ++tempBs)
               {
                    // First B-spline contribution
                    Real spanBSpline =
                        vec[(tempBs + kOrder - 1) - 1] - vec[tempBs - 1];
                    if (spanBSpline > NOD_THRESHOLD)
                    {
                         for (int inter = 0; inter <= kOrder - 2; ++inter)
                         {
                              Real spanPrior =
                                  vec[(tempBs + inter) - 1] - vec[tempBs - 1];
                              Real factorConst = spanPrior / spanBSpline;

                              // Constant term
                              for (int m = 0; m <= kOrder - 2; ++m)
                              {
                                   coe[idxCoe(m, inter, kOrder, tempBs)] +=
                                       coe[idxCoe(m, inter, kOrder - 1, tempBs)] *
                                       factorConst;
                              }

                              // Monomial term
                              Real spanCurrent =
                                  vec[(tempBs + inter + 1) - 1] -
                                  vec[(tempBs + inter) - 1];
                              Real factorMono = spanCurrent / spanBSpline;
                              for (int m = 1; m <= kOrder - 1; ++m)
                              {
                                   coe[idxCoe(m, inter, kOrder, tempBs)] +=
                                       coe[idxCoe(m - 1, inter, kOrder - 1, tempBs)] *
                                       factorMono;
                              }
                         }
                    }

                    // Second B-spline contribution
                    spanBSpline =
                        vec[(tempBs + kOrder) - 1] - vec[(tempBs + 1) - 1];
                    if (spanBSpline > NOD_THRESHOLD)
                    {
                         for (int inter = 1; inter <= kOrder - 1; ++inter)
                         {
                              Real spanToEnd =
                                  vec[(tempBs + kOrder) - 1] -
                                  vec[(tempBs + inter) - 1];
                              Real factorConst = spanToEnd / spanBSpline;

                              for (int m = 0; m <= kOrder - 2; ++m)
                              {
                                   coe[idxCoe(m, inter, kOrder, tempBs)] +=
                                       coe[idxCoe(m, inter - 1, kOrder - 1, tempBs + 1)] *
                                       factorConst;
                              }

                              Real spanCurrent =
                                  vec[(tempBs + inter + 1) - 1] -
                                  vec[(tempBs + inter) - 1];
                              Real factorMono = spanCurrent / spanBSpline;

                              for (int m = 1; m <= kOrder - 1; ++m)
                              {
                                   coe[idxCoe(m, inter, kOrder, tempBs)] -=
                                       coe[idxCoe(m - 1, inter - 1, kOrder - 1, tempBs + 1)] *
                                       factorMono;
                              }
                         }
                    }
               }
          }

          // Extract Coefficients = coe(:,:,order,1)
          for (int m = 0; m < order; ++m)
          {
               for (int inter = 0; inter < order; ++inter)
               {
                    coefficients[m + order * inter] =
                        coe[idxCoe(m, inter, order, 1)];
               }
          }
     }

     int BSpline::whichInterval(Real x) const
     {
          checkInitialization();

          int &i1 = saved_i1;
          int i2 = 0;

          if (x > gridAt(i1) && x <= gridAt(i1 + 1))
          {
               return i1;
          }
          if (x > gridAt(i1 + 1) &&
              x <= gridAt(std::min(i1 + 2, nNodes_ + order_ - 1)))
          {
               return i1 + 1;
          }

          i1 = 0;
          i2 = 0;

          if (x < gridAt(0))
          {
               return -1;
          }
          else if (x > gridAt(nNodes_ - 1))
          {
               return nNodes_;
          }

          if (x < gridAt(i1))
               i1 = 0;
          if (x > gridAt(i2))
               i2 = nNodes_ - 1;

          while (i2 - i1 > 1)
          {
               int mid = (i1 + i2) / 2;
               if (x > gridAt(mid))
               {
                    i1 = mid;
               }
               else
               {
                    i2 = mid;
               }
          }
          return i1;
     }

     Real BSpline::eval(Real x, int Bs, int derivativeOrder) const
     {
          checkInitialization();

          if (Bs < 1 || Bs > nBSplines_)
          {
               return 0.0;
          }

          int n = derivativeOrder;
          if (n < 0)
          {
               throw std::runtime_error("Negative derivative order in BSpline::eval");
          }
          if (n >= order_)
          {
               return 0.0;
          }

          int i = whichInterval(x);
          if (i < std::max(0, Bs - order_) ||
              i > std::min(nNodes_ - 1, Bs - 1))
          {
               return 0.0;
          }

          int j = i - Bs + order_;

          Real intervalWidth = gridAt(i + 1) - gridAt(i);
          Real oneOverInterval = 1.0 / intervalWidth;
          Real r = (x - gridAt(i)) * oneOverInterval;

          Real a = 1.0;
          Real w = 0.0;
          for (int k = n; k <= order_ - 1; ++k)
          {
               Real pf = partialFactorial(k, n);
               w += a * pf * coeff(k, j, Bs);
               a *= r;
          }

          Real oneOverIntervalToN = 1.0;
          for (int k = 1; k <= n; ++k)
          {
               oneOverIntervalToN *= oneOverInterval;
          }

          return w * oneOverIntervalToN;
     }

     Real BSpline::eval(Real x,
                        const Real *fc,
                        int fcSize,
                        int derivativeOrder,
                        bool skipFirst,
                        int Bsmin,
                        int Bsmax) const
     {
          checkInitialization();

          if (!fc || fcSize <= 0)
          {
               return 0.0;
          }

          int i = whichInterval(x);
          if (i < 0 || i > nNodes_ - 2)
          {
               return 0.0;
          }

          int Bsmi = 1;
          int nskip = 0;

          if (skipFirst)
          {
               Bsmi = 2;
               nskip = 1;
          }

          if (Bsmin >= 1)
          {
               Bsmi = std::max(Bsmin, Bsmi);
               nskip = std::max(Bsmi - 1, nskip);
          }

          Bsmi = std::max(Bsmi, i + 1);
          int Bsma = std::min(nBSplines_, i + order_);
          if (Bsmax >= 1)
          {
               Bsma = std::min(Bsmax, Bsma);
          }

          Real result = 0.0;
          for (int Bs = Bsmi; Bs <= Bsma; ++Bs)
          {
               int idxCoeff = Bs - nskip - 1; // 0-based
               if (idxCoeff < 0 || idxCoeff >= fcSize)
               {
                    continue;
               }
               Real c = fc[idxCoeff];
               if (c != 0.0)
               {
                    result += c * eval(x, Bs, derivativeOrder);
               }
          }

          return result;
     }

     Real BSpline::integral(const D2DFun &fun,
                            int Bs1,
                            int Bs2,
                            int braDerivativeOrder,
                            int ketDerivativeOrder,
                            const Real *parvec,
                            bool hasLowerBound,
                            Real lowerBound,
                            bool hasUpperBound,
                            Real upperBound,
                            bool hasBreakPoint,
                            Real breakPoint) const
     {
          checkInitialization();

          if (std::min(Bs1, Bs2) < 1 ||
              std::max(Bs1, Bs2) > nBSplines_ ||
              std::abs(Bs1 - Bs2) >= order_)
          {
               return 0.0;
          }

          Real a = gridAt(std::max(Bs1, Bs2) - order_);
          Real b = gridAt(std::min(Bs1, Bs2));

          if (hasLowerBound)
          {
               if (lowerBound >= b)
                    return 0.0;
               a = std::max(a, lowerBound);
          }
          if (hasUpperBound)
          {
               if (upperBound <= a)
                    return 0.0;
               b = std::min(b, upperBound);
          }

          int n1 = braDerivativeOrder;
          int n2 = ketDerivativeOrder;

          if (hasBreakPoint)
          {
               Real bp = breakPoint;
               const Real eps = std::numeric_limits<Real>::epsilon();
               if (bp <= a)
                    bp = a + eps;
               else if (bp >= b)
                    bp = b - eps;

               Real I1 = driverIntegral(fun, Bs1, Bs2, n1, n2, a, bp, parvec);
               Real I2 = driverIntegral(fun, Bs1, Bs2, n1, n2, bp, b, parvec);
               return I1 + I2;
          }
          else
          {
               return driverIntegral(fun, Bs1, Bs2, n1, n2, a, b, parvec);
          }
     }

     Real BSpline::driverIntegral(const D2DFun &fun,
                                  int Bs1, int Bs2,
                                  int n1, int n2,
                                  Real a, Real b,
                                  const Real *parvec) const
     {
          Real integral = 0.0;

          Real aPlus = upperLimitTo(a);
          Real bMinus = lowerLimitTo(b);

          int intervalMin = whichInterval(aPlus);
          int intervalMax = whichInterval(bMinus);

          for (int interval = intervalMin; interval <= intervalMax; ++interval)
          {
               Real lower = std::max(a, gridAt(interval));
               Real upper = std::min(b, gridAt(interval + 1));
               Real width = upper - lower;
               if (width < NOD_THRESHOLD)
                    continue;

               Real partialIntegral = 0.0;
               for (int iGauss = 0; iGauss < NGAUSS; ++iGauss)
               {
                    Real r = lower + width * gaussPoints[iGauss];
                    Real fVal = fun(r, parvec);
                    Real b1 = eval(r, Bs1, n1);
                    Real b2 = eval(r, Bs2, n2);
                    partialIntegral += fVal * b1 * b2 * gaussWeights[iGauss];
               }
               partialIntegral *= width;
               integral += partialIntegral;
          }

          return integral;
     }

} // namespace bspline
