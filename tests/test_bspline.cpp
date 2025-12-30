// tests/test_bspline.cpp
//
// Unit tests for the bspline::BSpline C++ translation of ModuleBSpline,
// using GoogleTest.
//
// These tests are meant to catch regressions as you evolve the code.
// They focus on:
//   - basic initialization and getters
//   - node-position mapping
//   - partition-of-unity property
//   - behavior outside of the knot interval
//   - consistency between "single-basis" evaluation and linear-combination
//   - integral symmetry and support overlap
//   - normalization factors

#include <gtest/gtest.h>

#include "BSpline.hpp"

#include <vector>
#include <cmath>
#include <algorithm>

// -----------------------------------------------------------------------------
// Helper local operators in the format expected by bspline::D2DFun
// -----------------------------------------------------------------------------

// Constant operator f(r) = 1
static double UnityOperator(double x, const double * /*parvec*/)
{
     (void)x;
     return 1.0;
}

// Linear operator f(r) = r
static double LinearOperator(double x, const double * /*parvec*/)
{
     return x;
}

// -----------------------------------------------------------------------------
// Test fixture: builds a small, uniform B-spline basis on [0,1]
// -----------------------------------------------------------------------------

class BSplineTest : public ::testing::Test
{
protected:
     using Real = bspline::Real;

     bspline::BSpline spline;
     int nNodes = 11; // small grid for tests
     int order = 4;   // cubic B-splines
     Real gridMin = 0.0;
     Real gridMax = 1.0;
     int nBSplines = 0;

     void SetUp() override
     {
          // Uniform grid on [gridMin, gridMax]
          std::vector<Real> grid(nNodes);
          for (int i = 0; i < nNodes; ++i)
          {
               grid[i] = gridMin + (gridMax - gridMin) * static_cast<Real>(i) / static_cast<Real>(nNodes - 1);
          }

          // Initialize B-spline basis
          int info = spline.init(nNodes, order, grid);
          ASSERT_EQ(info, 0) << "BSpline initialization failed";

          nBSplines = spline.getNBSplines();
          ASSERT_EQ(nBSplines, nNodes + order - 2)
              << "NBSplines should be NNodes + Order - 2 for standard knot multiplicities";
     }

     // Evaluate a *single* B-spline B_i(x) using the "linear combination" Eval:
     // set coeffs = e_i (unit vector) and call eval(x, coeffs.data(), nBSplines).
     Real evalSingle(Real x, int bsIndex) const
     {
          EXPECT_GE(bsIndex, 1);
          EXPECT_LE(bsIndex, nBSplines);

          if (bsIndex < 1 || bsIndex > nBSplines)
          {
               return Real(0); // avoid undefined behavior on bad input
          }

          std::vector<Real> coeffs(nBSplines, Real(0));
          coeffs[bsIndex - 1] = Real(1);

          return spline.eval(x, coeffs.data(), nBSplines);
     }
};

// -----------------------------------------------------------------------------
// Basic initialization and getters
// -----------------------------------------------------------------------------

TEST_F(BSplineTest, BasicInitializationAndGetters)
{
     EXPECT_EQ(spline.getNNodes(), nNodes);
     EXPECT_EQ(spline.getOrder(), order);
     EXPECT_EQ(spline.getNBSplines(), nBSplines);

     // Check that the "clamped" node positions match the uniform grid.
     // For a uniform [0,1] grid with nNodes=11:
     //   node 1   -> position 0.0
     //   node 2   -> position 0.1
     //   node 11  -> position 1.0
     bspline::Real pos1 = spline.getNodePosition(1);
     bspline::Real pos2 = spline.getNodePosition(2);
     bspline::Real pos11 = spline.getNodePosition(nNodes);

     EXPECT_NEAR(pos1, 0.0, 1e-14);
     EXPECT_NEAR(pos2, 0.1, 1e-14);
     EXPECT_NEAR(pos11, 1.0, 1e-14);
}

// -----------------------------------------------------------------------------
// Partition-of-unity: sum_i B_i(x) ≈ 1 in the interior of the knot span
// -----------------------------------------------------------------------------

TEST_F(BSplineTest, PartitionOfUnity)
{
     // Sample a few interior points. We avoid exactly the endpoints to
     // stay away from potential edge degeneracies.
     std::vector<Real> xs = {0.05, 0.25, 0.5, 0.75, 0.95};

     for (Real x : xs)
     {
          Real sum = 0.0;
          for (int i = 1; i <= nBSplines; ++i)
          {
               sum += evalSingle(x, i);
          }
          EXPECT_NEAR(sum, 1.0, 1e-11) << "Partition of unity violated at x=" << x;
     }
}

// -----------------------------------------------------------------------------
// Outside the knot interval, all B-splines should evaluate to ~0
// -----------------------------------------------------------------------------

TEST_F(BSplineTest, ValuesZeroOutsideSupport)
{
     const Real xLeft = gridMin - 0.1;  // a bit to the left of the domain
     const Real xRight = gridMax + 0.1; // a bit to the right of the domain

     for (int i = 1; i <= nBSplines; ++i)
     {
          Real vLeft = evalSingle(xLeft, i);
          Real vRight = evalSingle(xRight, i);

          EXPECT_NEAR(vLeft, 0.0, 1e-14) << "B-spline " << i << " not ~0 at left of domain";
          EXPECT_NEAR(vRight, 0.0, 1e-14) << "B-spline " << i << " not ~0 at right of domain";
     }
}

// -----------------------------------------------------------------------------
// Normalization factors: f(i)^2 * ∫ B_i(r)^2 r^2 dr ≈ 1
// -----------------------------------------------------------------------------

TEST_F(BSplineTest, NormFactorsProduceUnitL2Norm)
{
     for (int i = 1; i <= nBSplines; ++i)
     {
          double Iii = spline.integral(&UnityOperator, i, i); // ∫ B_i B_i * 1 * r^2 dr
          double f = spline.getNormFactor(i);
          double normSq = f * f * Iii;

          EXPECT_NEAR(normSq, 1.0, 1e-11)
              << "Normalized B-spline " << i << " does not have unit L2 norm";
     }
}

// -----------------------------------------------------------------------------
// Function evaluation with coefficient vector should match explicit sum:
//   f(x) = ∑_i c_i B_i(x)
// -----------------------------------------------------------------------------

TEST_F(BSplineTest, FunctionEvalMatchesExplicitLinearCombination)
{
     using Real = bspline::Real;

     // Set some simple, deterministic coefficients
     std::vector<Real> coeffs(nBSplines);
     for (int i = 0; i < nBSplines; ++i)
     {
          coeffs[i] = static_cast<Real>(i + 1) * 0.1; // 0.1, 0.2, 0.3, ...
     }

     std::vector<Real> xs = {0.1, 0.33, 0.57, 0.9};

     for (Real x : xs)
     {
          // Evaluation via vector interface
          Real f_vec = spline.eval(x, coeffs.data(), nBSplines);

          // Explicit sum of c_i * B_i(x)
          Real f_sum = 0.0;
          for (int i = 1; i <= nBSplines; ++i)
          {
               f_sum += coeffs[i - 1] * evalSingle(x, i);
          }

          EXPECT_NEAR(f_vec, f_sum, 1e-11)
              << "Function eval mismatch at x=" << x;
     }
}

// -----------------------------------------------------------------------------
// Integral symmetry: for real B-splines and real local operators,
// ∫ B_i f B_j r^2 dr should be symmetric in (i,j).
// -----------------------------------------------------------------------------

TEST_F(BSplineTest, IntegralSymmetryWithUnity)
{
     for (int i = 1; i <= nBSplines; ++i)
     {
          for (int j = i; j <= std::min(i + order - 1, nBSplines); ++j)
          {
               double Iij = spline.integral(&UnityOperator, i, j);
               double Iji = spline.integral(&UnityOperator, j, i);

               EXPECT_NEAR(Iij, Iji, 1e-11)
                   << "Integral not symmetric for indices (" << i << "," << j << ")";
          }
     }
}

// Same symmetry test but with a simple non-constant local operator f(r) = r.
TEST_F(BSplineTest, IntegralSymmetryWithLinearWeight)
{
     for (int i = 1; i <= nBSplines; ++i)
     {
          for (int j = i; j <= std::min(i + order - 1, nBSplines); ++j)
          {
               double Iij = spline.integral(&LinearOperator, i, j);
               double Iji = spline.integral(&LinearOperator, j, i);

               EXPECT_NEAR(Iij, Iji, 1e-11)
                   << "Integral with linear weight not symmetric for indices ("
                   << i << "," << j << ")";
          }
     }
}

// -----------------------------------------------------------------------------
// Support overlap: if |i - j| >= order, supports do not overlap and
// the integral must be (numerically) zero.
// -----------------------------------------------------------------------------

TEST_F(BSplineTest, IntegralZeroForNonOverlappingBsplines)
{
     for (int i = 1; i <= nBSplines; ++i)
     {
          int j = i + order; // this should be outside the overlap band
          if (j > nBSplines)
               break;

          double Iij = spline.integral(&UnityOperator, i, j);
          EXPECT_NEAR(Iij, 0.0, 1e-12)
              << "Integral between non-overlapping B-splines ("
              << i << "," << j << ") is not ~0";
     }
}

// -----------------------------------------------------------------------------
// (Optional) Quick sanity: integrals are positive on the diagonal with Unity
// -----------------------------------------------------------------------------

TEST_F(BSplineTest, DiagonalIntegralsArePositive)
{
     for (int i = 1; i <= nBSplines; ++i)
     {
          double Iii = spline.integral(&UnityOperator, i, i);
          EXPECT_GT(Iii, 0.0) << "Self integral of B-spline " << i << " is not positive";
     }
}
