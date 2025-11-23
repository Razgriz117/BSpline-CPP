#pragma once

#include <vector>
#include <functional>

// B-spline implementation in C++ that mirrors the behavior of
// Luca Argenti's Fortran ModuleBSpline for the use-case in Template.f90.
namespace bspline
{

     // Scalar type used throughout (equivalent to double precision in Fortran).
     using Real = double;

     // Type of local operator functions f(x, parvec):
     //   - x:     radial coordinate
     //   - parvec: optional pointer to extra parameters (can be nullptr)
     using D2DFun = std::function<Real(Real x, const Real *parvec)>;

     // Forward declaration of the internal unity function used for normalization.
     // This is logically the same as the module's Unity in the Fortran code,
     // but distinct from any Unity defined in user code (e.g. main.cpp).
     Real unityFunction(Real x, const Real *parvec);

     // C++ version of ClassBSpline from ModuleBSpline.
     class BSpline
     {
     public:
          BSpline() = default;
          ~BSpline() = default;

          // Disable copy assignment/constructor only if you want to manage explicitly;
          // here we allow default copying (std::vector has deep copy semantics),
          // which mirrors Fortran's ASSIGNMENT(=) copy of ClassBSpline.
          BSpline(const BSpline &) = default;
          BSpline &operator=(const BSpline &) = default;

          // Initialize a B-spline basis set.
          //
          // Arguments:
          //   numberOfNodes : number of *distinct* grid nodes (>= 2)
          //   order         : spline order k (>= 1)
          //   grid          : vector of grid points (size >= numberOfNodes),
          //                   non-decreasing, representing the distinct nodes.
          //
          // Return value (same semantics as Fortran's IOSTAT):
          //   0   : success
          //  -1   : invalid numberOfNodes (< 2)
          //  -2   : invalid order (< 1)
          //  -4   : grid.size() < numberOfNodes
          //   1   : grid not non-decreasing
          int init(int numberOfNodes, int order, const std::vector<Real> &grid);

          // Free allocated data (equivalent to BSplineFree in Fortran).
          void free();

          // Accessors (equivalent to GetNNodes, GetOrder, GetNBSplines).
          int getNNodes() const { return nNodes_; }
          int getOrder() const { return order_; }
          int getNBSplines() const { return nBSplines_; }

          // Evaluate the n-th derivative (n = 0, 1, ...) of a single B-spline Bs at x.
          //
          //  - Bs is a 1-based index: 1 <= Bs <= nBSplines_ (mirrors Fortran)
          //  - derivativeOrder is the order of derivative (n in the Fortran code)
          //
          // Returns 0 if Bs is out of range or derivativeOrder >= order_ (as in Fortran).
          Real eval(Real x, int Bs, int derivativeOrder = 0) const;

          // Evaluate the k-th derivative of a function expressed as a linear combination
          // of B-splines:
          //
          //   f(x) = sum_{Bs = Bsmin..Bsmax} fc(Bs - nskip) * B_spline_Bs(x)
          //
          // This mirrors BSplineFunctionEval in the Fortran module.
          //
          // Arguments:
          //   x               : point at which to evaluate
          //   fc              : pointer to coefficient array (Fortran: fc(:))
          //   fcSize          : number of elements in fc
          //   derivativeOrder : derivative order k (default 0)
          //   skipFirst       : if true, coefficients start from Bs=2 (Fortran SKIP_FIRST)
          //   Bsmin, Bsmax    : optional bounds on B-spline indices to include
          //
          // For the use in Template.f90, we call this with:
          //   - skipFirst = false
          //   - Bsmin, Bsmax = -1 (meaning "not provided")
          Real eval(Real x,
                    const Real *fc,
                    int fcSize,
                    int derivativeOrder = 0,
                    bool skipFirst = false,
                    int Bsmin = -1,
                    int Bsmax = -1) const;

          // Compute the integral
          //
          //   ∫_a^b [ d^{n1} Bs_i(r)/dr^{n1} ] * f(r) *
          //         [ d^{n2} Bs_j(r)/dr^{n2} ] * r^2 dr
          //
          // where f(r) is a local operator.
          //
          // This is a C++ analog of BSplineIntegral in Fortran, but simplified
          // to cover the use-cases in Template.f90:
          //   - optional derivative orders n1, n2 (default 0)
          //   - optional parvec (can be nullptr)
          //   - optional lower / upper bounds (flags indicate whether they're set)
          //   - optional breakpoint (for splitting, not used in Template.f90)
          //
          // For Template.f90 we only need:
          //   integral(f, Bs1, Bs2)                     -- overlap
          //   integral(f, Bs1, Bs2, 1, 1)               -- kinetic part
          //   integral(f, Bs1, Bs2, 0, 0, parvec)       -- potential part
          Real integral(const D2DFun &fun,
                        int Bs1,
                        int Bs2,
                        int braDerivativeOrder = 0,
                        int ketDerivativeOrder = 0,
                        const Real *parvec = nullptr,
                        bool hasLowerBound = false,
                        Real lowerBound = 0.0,
                        bool hasUpperBound = false,
                        Real upperBound = 0.0,
                        bool hasBreakPoint = false,
                        Real breakPoint = 0.0) const;

     private:
          // Internal state (mirrors Fortran type components).
          bool initialized_ = false; // Equivalent of s%INITIALIZED
          int nNodes_ = 0;           // Number of distinct nodes
          int order_ = 0;            // Spline order
          int nBSplines_ = 0;        // Total number of B-splines = nNodes_ + order_ - 2

          // Extended grid: corresponds to s%Grid(-(order-1) : nNodes_-1 + order_-1)
          //
          // This holds:
          //   - repeated first node order_-1 times at the beginning
          //   - original nodes
          //   - repeated last node order_-1 times at the end
          //
          // We emulate Fortran's negative/shifted indices using an offset.
          std::vector<Real> grid_; // size = nNodes_ + 2*order_ - 2

          // Normalization factors f(i) = 1 / ||B_i||_L2
          std::vector<Real> normalization_; // size = nBSplines_

          // Polynomial coefficients c(m, j, Bs) of B-splines.
          //
          // Fortran dimensions:
          //   c(0:order_-1, 0:order_-1, 1:nBSplines_)
          //
          // We store them flattened as:
          //   index = m + order_ * (j + order_ * (Bs-1))
          // where:
          //   m   = 0..order_-1   (polynomial degree)
          //   j   = 0..order_-1   (interval index in support)
          //   Bs  = 1..nBSplines_ (B-spline index)
          std::vector<Real> coeffs_;

          // Helper for mapping a "Fortran-style" grid index (which can be negative)
          // to the actual vector index in grid_.
          inline int idxGrid(int fortranIndex) const
          {
               // In Fortran: lower bound is -(order_-1); we shift by + (order_-1)
               return fortranIndex + (order_ - 1);
          }

          inline Real gridAt(int fortranIndex) const
          {
               return grid_[idxGrid(fortranIndex)];
          }

          inline Real &gridAt(int fortranIndex)
          {
               return grid_[idxGrid(fortranIndex)];
          }

          // Accessor for coefficient c(m, j, Bs).
          inline Real coeff(int m, int j, int Bs) const
          {
               const int idx = m + order_ * (j + order_ * (Bs - 1));
               return coeffs_[idx];
          }

          inline Real &coeff(int m, int j, int Bs)
          {
               const int idx = m + order_ * (j + order_ * (Bs - 1));
               return coeffs_[idx];
          }

          // Compute coefficients of a single B-spline given order and a vector of
          // (order+1) consecutive nodes vec(1:order+1). This mirrors
          // ComputeCoefficientsSingleBSpline in Fortran.
          //
          // Arguments:
          //   order        : spline order
          //   vec          : pointer to an array of size order+1
          //   coefficients : output array of size order * order
          void computeCoefficientsSingleBSpline(int order,
                                                const Real *vec,
                                                Real *coefficients) const;

          // Find which interval (g_i, g_{i+1}] contains x.
          // This mirrors which_interval(x, s) in Fortran.
          //
          // Returns:
          //   i in [0, nNodes_-1] such that x ∈ (g_i, g_{i+1}]
          //   -1 if x < g_0
          //   nNodes_ if x > g_{nNodes_-1}
          int whichInterval(Real x) const;

          // Internal driver for evaluating integrals once parameters are checked.
          // This corresponds to BSplineDriverIntegral in Fortran.
          Real driverIntegral(const D2DFun &fun,
                              int Bs1, int Bs2,
                              int n1, int n2,
                              Real a, Real b,
                              const Real *parvec) const;

          // Ensure that init() has been called before using this object.
          void checkInitialization() const;
     };

} // namespace bspline
