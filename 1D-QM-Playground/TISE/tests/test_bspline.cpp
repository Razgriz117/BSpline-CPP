#include <gtest/gtest.h>
#include "BSpline.hpp"

#include <algorithm>
#include <vector>

// ---------------------------------------------------------------------------
// Local operator helpers (D2DFun signature)
// ---------------------------------------------------------------------------

static double UnityOp(double /*x*/, const double * /*p*/) { return 1.0; }
static double LinearOp(double x, const double * /*p*/) { return x; }
static double ScaledUnityOp(double /*x*/, const double *p) { return p[0]; }

// ---------------------------------------------------------------------------
// Grid builders
// ---------------------------------------------------------------------------

static std::vector<double> uniformGrid(int nNodes, double lo, double hi)
{
    std::vector<double> g(nNodes);
    for (int i = 0; i < nNodes; ++i)
        g[i] = lo + (hi - lo) * static_cast<double>(i) / (nNodes - 1);
    return g;
}

static std::vector<double> quadraticGrid(int nNodes, double lo, double hi)
{
    std::vector<double> g(nNodes);
    for (int i = 0; i < nNodes; ++i)
    {
        double t = static_cast<double>(i) / (nNodes - 1);
        g[i] = lo + (hi - lo) * t * t;
    }
    return g;
}

// ---------------------------------------------------------------------------
// Primary test fixture: uniform grid on [0,1], nNodes=11, order=4
// ---------------------------------------------------------------------------

class BSplineTest : public ::testing::Test
{
protected:
    using Real = bspline::Real;

    bspline::BSpline spline;
    int nNodes    = 11;
    int order     = 4;
    Real gridMin  = 0.0;
    Real gridMax  = 1.0;
    int nBSplines = 0;

    void SetUp() override
    {
        auto grid = uniformGrid(nNodes, gridMin, gridMax);
        int info  = spline.init(nNodes, order, grid);
        ASSERT_EQ(info, 0) << "BSpline initialization failed";
        nBSplines = spline.getNBSplines();
        ASSERT_EQ(nBSplines, nNodes + order - 2);
    }

    Real evalSingle(Real x, int bsIndex) const
    {
        return spline.eval(x, bsIndex);
    }

    Real evalViaCoeffs(Real x, int bsIndex) const
    {
        std::vector<Real> c(nBSplines, 0.0);
        c[bsIndex - 1] = 1.0;
        return spline.eval(x, c.data(), nBSplines);
    }
};

// ---------------------------------------------------------------------------
// Section 1: init() return codes
// ---------------------------------------------------------------------------

TEST(BSplineInitTest, InvalidNNodes_ReturnsMinus1)
{
    bspline::BSpline s;
    auto grid = uniformGrid(2, 0.0, 1.0);
    EXPECT_EQ(s.init(1, 4, grid), -1);
    EXPECT_EQ(s.init(0, 4, grid), -1);
}

TEST(BSplineInitTest, InvalidOrder_ReturnsMinus2)
{
    bspline::BSpline s;
    auto grid = uniformGrid(5, 0.0, 1.0);
    EXPECT_EQ(s.init(5,  0, grid), -2);
    EXPECT_EQ(s.init(5, -1, grid), -2);
}

TEST(BSplineInitTest, GridTooSmall_ReturnsMinus4)
{
    bspline::BSpline s;
    std::vector<double> grid = {0.0, 0.5, 1.0}; // size=3, nNodes=5
    EXPECT_EQ(s.init(5, 4, grid), -4);
}

TEST(BSplineInitTest, NonMonotonicGrid_ReturnsOne)
{
    bspline::BSpline s;
    std::vector<double> grid = {0.0, 0.5, 0.3, 0.8, 1.0};
    EXPECT_EQ(s.init(5, 4, grid), 1);
}

TEST(BSplineInitTest, ValidInit_ReturnsZero)
{
    bspline::BSpline s;
    auto grid = uniformGrid(11, 0.0, 1.0);
    EXPECT_EQ(s.init(11, 4, grid), 0);
}

// ---------------------------------------------------------------------------
// Section 2: free() resets state
// ---------------------------------------------------------------------------

TEST(BSplineFreeTest, FreeResetsGetters)
{
    bspline::BSpline s;
    auto grid = uniformGrid(11, 0.0, 1.0);
    ASSERT_EQ(s.init(11, 4, grid), 0);
    EXPECT_NE(s.getNBSplines(), 0);

    s.free();

    EXPECT_EQ(s.getNNodes(),    0);
    EXPECT_EQ(s.getOrder(),     0);
    EXPECT_EQ(s.getNBSplines(), 0);
}

TEST(BSplineFreeTest, ReinitAfterFreeWorks)
{
    bspline::BSpline s;
    auto g1 = uniformGrid(11, 0.0, 1.0);
    ASSERT_EQ(s.init(11, 4, g1), 0);
    s.free();

    auto g2 = uniformGrid(7, 0.0, 2.0);
    EXPECT_EQ(s.init(7, 6, g2), 0);
    EXPECT_EQ(s.getNBSplines(), 7 + 6 - 2);
}

// ---------------------------------------------------------------------------
// Section 3: Basic getters
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, BasicInitializationAndGetters)
{
    EXPECT_EQ(spline.getNNodes(),    nNodes);
    EXPECT_EQ(spline.getOrder(),     order);
    EXPECT_EQ(spline.getNBSplines(), nBSplines);
    EXPECT_EQ(nBSplines, nNodes + order - 2);
}

// ---------------------------------------------------------------------------
// Section 4: eval(x, Bs) — out-of-range Bs
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, EvalOutOfRangeBsReturnsZero)
{
    EXPECT_EQ(spline.eval(0.5, 0),             0.0);
    EXPECT_EQ(spline.eval(0.5, -1),            0.0);
    EXPECT_EQ(spline.eval(0.5, nBSplines + 1), 0.0);
}

// ---------------------------------------------------------------------------
// Section 5: Partition of unity
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, PartitionOfUnity)
{
    std::vector<double> xs = {0.05, 0.25, 0.5, 0.75, 0.95};
    for (double x : xs)
    {
        double sum = 0.0;
        for (int i = 1; i <= nBSplines; ++i)
            sum += evalSingle(x, i);
        EXPECT_NEAR(sum, 1.0, 1e-11)
            << "Partition of unity violated at x=" << x;
    }
}

// ---------------------------------------------------------------------------
// Section 6: Derivative of partition of unity = 0
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, DerivativePartitionOfUnityIsZero)
{
    std::vector<double> xs = {0.1, 0.3, 0.5, 0.7, 0.9};
    for (double x : xs)
    {
        double dsum = 0.0;
        for (int i = 1; i <= nBSplines; ++i)
            dsum += spline.eval(x, i, /*derivativeOrder=*/1);
        EXPECT_NEAR(dsum, 0.0, 1e-10)
            << "Derivative of partition of unity != 0 at x=" << x;
    }
}

// ---------------------------------------------------------------------------
// Section 7: Values zero outside domain
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, ValuesZeroOutsideSupport)
{
    const double xLeft  = gridMin - 0.1;
    const double xRight = gridMax + 0.1;
    for (int i = 1; i <= nBSplines; ++i)
    {
        EXPECT_NEAR(evalSingle(xLeft,  i), 0.0, 1e-14)
            << "B-spline " << i << " nonzero left of domain";
        EXPECT_NEAR(evalSingle(xRight, i), 0.0, 1e-14)
            << "B-spline " << i << " nonzero right of domain";
    }
}

// ---------------------------------------------------------------------------
// Section 8: Linear combination eval
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, FunctionEvalMatchesExplicitLinearCombination)
{
    std::vector<double> coeffs(nBSplines);
    for (int i = 0; i < nBSplines; ++i)
        coeffs[i] = (i + 1) * 0.1;

    std::vector<double> xs = {0.1, 0.33, 0.57, 0.9};
    for (double x : xs)
    {
        double f_vec = spline.eval(x, coeffs.data(), nBSplines);
        double f_sum = 0.0;
        for (int i = 1; i <= nBSplines; ++i)
            f_sum += coeffs[i - 1] * evalSingle(x, i);
        EXPECT_NEAR(f_vec, f_sum, 1e-11)
            << "Function eval mismatch at x=" << x;
    }
}

TEST_F(BSplineTest, SingleBasisEvalAndCoeffEvalAgree)
{
    std::vector<double> xs = {0.05, 0.4, 0.6, 0.95};
    for (double x : xs)
    {
        for (int i = 1; i <= nBSplines; ++i)
        {
            EXPECT_NEAR(evalSingle(x, i), evalViaCoeffs(x, i), 1e-14)
                << "Mismatch for Bs=" << i << " at x=" << x;
        }
    }
}

// ---------------------------------------------------------------------------
// Section 9: eval with skipFirst
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, SkipFirstExcludesFirstBSpline)
{
    std::vector<double> ones(nBSplines, 1.0);
    std::vector<double> xs = {0.2, 0.5, 0.8};
    for (double x : xs)
    {
        double withSkip = spline.eval(x, ones.data(), nBSplines,
                                      /*derivOrder=*/0, /*skipFirst=*/true);
        double explicit_sum = 0.0;
        for (int i = 2; i <= nBSplines; ++i)
            explicit_sum += evalSingle(x, i);
        EXPECT_NEAR(withSkip, explicit_sum, 1e-11)
            << "skipFirst mismatch at x=" << x;
    }
}

// ---------------------------------------------------------------------------
// Section 10: eval with Bsmin / Bsmax bounds
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, BsminBsmaxLimitsSummation)
{
    std::vector<double> ones(nBSplines, 1.0);
    double x = 0.5;

    double with_bounds = spline.eval(x, ones.data(), nBSplines,
                                     /*derivOrder=*/0,
                                     /*skipFirst=*/false,
                                     /*Bsmin=*/3,
                                     /*Bsmax=*/5);
    double explicit_sum = 0.0;
    for (int i = 3; i <= 5; ++i)
        explicit_sum += evalSingle(x, i);

    EXPECT_NEAR(with_bounds, explicit_sum, 1e-11);
}

// ---------------------------------------------------------------------------
// Section 11: Self-overlap integral (positive normalization sanity)
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, SelfOverlapIntegralIsPositive)
{
    for (int i = 1; i <= nBSplines; ++i)
    {
        double I = spline.integral(&UnityOp, i, i);
        EXPECT_GT(I, 0.0)
            << "Self-overlap integral of B-spline " << i << " is not positive";
    }
}

// ---------------------------------------------------------------------------
// Section 12: Integral symmetry
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, IntegralSymmetryWithUnity)
{
    for (int i = 1; i <= nBSplines; ++i)
    {
        for (int j = i; j <= std::min(i + order - 1, nBSplines); ++j)
        {
            double Iij = spline.integral(&UnityOp, i, j);
            double Iji = spline.integral(&UnityOp, j, i);
            EXPECT_NEAR(Iij, Iji, 1e-11)
                << "Unity integral not symmetric for (" << i << "," << j << ")";
        }
    }
}

TEST_F(BSplineTest, IntegralSymmetryWithLinearWeight)
{
    for (int i = 1; i <= nBSplines; ++i)
    {
        for (int j = i; j <= std::min(i + order - 1, nBSplines); ++j)
        {
            double Iij = spline.integral(&LinearOp, i, j);
            double Iji = spline.integral(&LinearOp, j, i);
            EXPECT_NEAR(Iij, Iji, 1e-11)
                << "Linear integral not symmetric for (" << i << "," << j << ")";
        }
    }
}

// ---------------------------------------------------------------------------
// Section 13: Zero integral for non-overlapping B-splines
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, IntegralZeroForNonOverlappingBSplines)
{
    for (int i = 1; i <= nBSplines; ++i)
    {
        int j = i + order;
        if (j > nBSplines) break;
        double Iij = spline.integral(&UnityOp, i, j);
        EXPECT_NEAR(Iij, 0.0, 1e-12)
            << "Non-overlapping integral (" << i << "," << j << ") is not ~0";
    }
}

// ---------------------------------------------------------------------------
// Section 14: Kinetic integrals (derivative orders 1,1)
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, KineticIntegralSymmetry)
{
    for (int i = 1; i <= nBSplines; ++i)
    {
        for (int j = i; j <= std::min(i + order - 1, nBSplines); ++j)
        {
            double Tij = spline.integral(&UnityOp, i, j, 1, 1);
            double Tji = spline.integral(&UnityOp, j, i, 1, 1);
            EXPECT_NEAR(Tij, Tji, 1e-11)
                << "Kinetic integral not symmetric for (" << i << "," << j << ")";
        }
    }
}

TEST_F(BSplineTest, KineticDiagonalIsPositive)
{
    for (int i = 1; i <= nBSplines; ++i)
    {
        double T = spline.integral(&UnityOp, i, i, 1, 1);
        EXPECT_GT(T, 0.0)
            << "Diagonal kinetic integral for B-spline " << i << " is not positive";
    }
}

// ---------------------------------------------------------------------------
// Section 15: integral() with parvec
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, IntegralWithParvecScalesResult)
{
    double scale = 3.14;
    const double *pv = &scale;

    for (int i = 1; i <= nBSplines; ++i)
    {
        double I_unity  = spline.integral(&UnityOp,       i, i);
        double I_scaled = spline.integral(&ScaledUnityOp, i, i, 0, 0, pv);
        EXPECT_NEAR(I_scaled, scale * I_unity, 1e-10)
            << "Parvec-scaled integral mismatch for B-spline " << i;
    }
}

// ---------------------------------------------------------------------------
// Section 16: integral() with lower/upper bound constraints
// ---------------------------------------------------------------------------

TEST_F(BSplineTest, IntegralWithUpperBound_SmallerThanFull)
{
    for (int i = 2; i <= nBSplines - 1; ++i)
    {
        double I_full = spline.integral(&UnityOp, i, i);
        double I_half = spline.integral(&UnityOp, i, i,
                                        0, 0, nullptr,
                                        /*hasLowerBound=*/false, 0.0,
                                        /*hasUpperBound=*/true, 0.5);
        EXPECT_LE(I_half, I_full + 1e-14)
            << "Upper-bound integral exceeds full integral for B-spline " << i;
        EXPECT_GE(I_half, 0.0);
    }
}

TEST_F(BSplineTest, IntegralLowerPlusUpperEqualsFullForMiddleBSpline)
{
    int i = nBSplines / 2;
    double I_full  = spline.integral(&UnityOp, i, i);
    double I_left  = spline.integral(&UnityOp, i, i,
                                     0, 0, nullptr,
                                     false, 0.0, true, 0.5);
    double I_right = spline.integral(&UnityOp, i, i,
                                     0, 0, nullptr,
                                     true, 0.5, false, 0.0);
    EXPECT_NEAR(I_left + I_right, I_full, 1e-11)
        << "Split integral does not sum to full for B-spline " << i;
}

// ---------------------------------------------------------------------------
// Section 17: Multiple orders (parameterized)
// ---------------------------------------------------------------------------

class BSplineOrderTest : public ::testing::TestWithParam<int> {};

TEST_P(BSplineOrderTest, PartitionOfUnityHoldsForOrder)
{
    int ord = GetParam();
    int nN  = std::max(4, ord + 2);
    bspline::BSpline s;
    auto g = uniformGrid(nN, 0.0, 1.0);
    ASSERT_EQ(s.init(nN, ord, g), 0) << "Init failed for order=" << ord;

    int nBs = s.getNBSplines();
    ASSERT_EQ(nBs, nN + ord - 2);

    std::vector<double> xs = {0.1, 0.4, 0.6, 0.9};
    for (double x : xs)
    {
        double sum = 0.0;
        for (int i = 1; i <= nBs; ++i)
            sum += s.eval(x, i);
        EXPECT_NEAR(sum, 1.0, 1e-10)
            << "Partition of unity failed for order=" << ord << " at x=" << x;
    }
}

INSTANTIATE_TEST_SUITE_P(Orders, BSplineOrderTest,
    ::testing::Values(1, 4, 8, 12));

// ---------------------------------------------------------------------------
// Section 18: Non-uniform (quadratic) grid
// ---------------------------------------------------------------------------

TEST(BSplineNonUniformGridTest, PartitionOfUnityOnQuadraticGrid)
{
    bspline::BSpline s;
    auto g = quadraticGrid(15, 0.0, 1.0);
    ASSERT_EQ(s.init(15, 4, g), 0);

    int nBs = s.getNBSplines();
    std::vector<double> xs = {0.02, 0.1, 0.3, 0.7, 0.95};
    for (double x : xs)
    {
        double sum = 0.0;
        for (int i = 1; i <= nBs; ++i)
            sum += s.eval(x, i);
        EXPECT_NEAR(sum, 1.0, 1e-11)
            << "Partition of unity on non-uniform grid failed at x=" << x;
    }
}

TEST(BSplineNonUniformGridTest, IntegralSymmetryOnNonUniformGrid)
{
    bspline::BSpline s;
    auto g = quadraticGrid(15, 0.0, 1.0);
    ASSERT_EQ(s.init(15, 4, g), 0);

    int nBs   = s.getNBSplines();
    int order = s.getOrder();
    for (int i = 1; i <= nBs; ++i)
    {
        for (int j = i; j <= std::min(i + order - 1, nBs); ++j)
        {
            double Iij = s.integral(&UnityOp, i, j);
            double Iji = s.integral(&UnityOp, j, i);
            EXPECT_NEAR(Iij, Iji, 1e-11)
                << "Non-uniform grid integral not symmetric for ("
                << i << "," << j << ")";
        }
    }
}
