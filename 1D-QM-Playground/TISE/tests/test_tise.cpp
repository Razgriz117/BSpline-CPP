#include <gtest/gtest.h>
#include "tise.hpp"
#include "BSpline.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

// ---------------------------------------------------------------------------
// buildUniformRadialGrid
// ---------------------------------------------------------------------------

TEST(BuildUniformRadialGridTest, SizeIsNNodes)
{
    auto g = tise::buildUniformRadialGrid(11, 0.0, 1.0);
    EXPECT_EQ(static_cast<int>(g.size()), 11);
}

TEST(BuildUniformRadialGridTest, FirstAndLastValues)
{
    auto g = tise::buildUniformRadialGrid(11, 0.0, 10.0);
    EXPECT_NEAR(g.front(), 0.0,  1e-15);
    EXPECT_NEAR(g.back(),  10.0, 1e-15);
}

TEST(BuildUniformRadialGridTest, UniformSpacing)
{
    int n = 6;
    auto g = tise::buildUniformRadialGrid(n, 0.0, 5.0);
    double expected_dx = 1.0; // 5.0 / (6-1)
    for (int i = 1; i < n; ++i)
        EXPECT_NEAR(g[i] - g[i-1], expected_dx, 1e-14);
}

TEST(BuildUniformRadialGridTest, TwoNodeGrid)
{
    auto g = tise::buildUniformRadialGrid(2, 3.0, 7.0);
    ASSERT_EQ(static_cast<int>(g.size()), 2);
    EXPECT_NEAR(g[0], 3.0, 1e-15);
    EXPECT_NEAR(g[1], 7.0, 1e-15);
}

// ---------------------------------------------------------------------------
// radialPotential
// ---------------------------------------------------------------------------

TEST(RadialPotentialTest, L0IsMinusOneOverX)
{
    // L=0: V(x) = -1/x
    EXPECT_NEAR(tise::radialPotential(1.0,  0), -1.0,   1e-15);
    EXPECT_NEAR(tise::radialPotential(2.0,  0), -0.5,   1e-15);
    EXPECT_NEAR(tise::radialPotential(0.5,  0), -2.0,   1e-15);
}

TEST(RadialPotentialTest, L1KnownValues)
{
    // L=1: V(x) = l*(l+1)/(2x^2) - 1/x = 1*2/(2x^2) - 1/x = 1/x^2 - 1/x
    // at x=1: 1 - 1 = 0
    EXPECT_NEAR(tise::radialPotential(1.0, 1),  0.0,   1e-15);
    // at x=2: 1/4 - 1/2 = -0.25
    EXPECT_NEAR(tise::radialPotential(2.0, 1), -0.25,  1e-15);
}

TEST(RadialPotentialTest, L2KnownValues)
{
    // L=2: V(x) = 3/x^2 - 1/x
    // at x=1: 3 - 1 = 2
    EXPECT_NEAR(tise::radialPotential(1.0, 2), 2.0,     1e-15);
}

// ---------------------------------------------------------------------------
// analyticHydrogenEnergy
// ---------------------------------------------------------------------------

TEST(AnalyticHydrogenEnergyTest, GroundStateL0)
{
    // n=1, L=0: E = -1/(2*1^2) = -0.5
    EXPECT_NEAR(tise::analyticHydrogenEnergy(1, 0), -0.5, 1e-15);
}

TEST(AnalyticHydrogenEnergyTest, GroundStateL1)
{
    // n=1, L=1: n_eff=2, E = -1/(2*4) = -0.125
    EXPECT_NEAR(tise::analyticHydrogenEnergy(1, 1), -0.125, 1e-15);
}

TEST(AnalyticHydrogenEnergyTest, ExcitedStates)
{
    // n=2, L=0: n_eff=2, E = -1/8 = -0.125
    EXPECT_NEAR(tise::analyticHydrogenEnergy(2, 0), -0.125, 1e-15);
    // n=3, L=0: n_eff=3, E = -1/18
    EXPECT_NEAR(tise::analyticHydrogenEnergy(3, 0), -1.0/18.0, 1e-15);
}

TEST(AnalyticHydrogenEnergyTest, IncreaseWithN)
{
    // Energies are negative and approach 0 from below as n increases,
    // so E(n+1) > E(n) (less negative = higher energy).
    for (int n = 1; n < 5; ++n)
        EXPECT_GT(tise::analyticHydrogenEnergy(n+1, 0),
                  tise::analyticHydrogenEnergy(n,   0));
}

// ---------------------------------------------------------------------------
// eigenvalueError
// ---------------------------------------------------------------------------

TEST(EigenvalueErrorTest, ExactValueGivesZero)
{
    double exact = tise::analyticHydrogenEnergy(1, 0);
    EXPECT_NEAR(tise::eigenvalueError(exact, 1, 0), 0.0, 1e-15);
}

TEST(EigenvalueErrorTest, PositiveForEnergyAboveExact)
{
    double exact = tise::analyticHydrogenEnergy(1, 0);
    EXPECT_GT(tise::eigenvalueError(exact + 0.001, 1, 0), 0.0);
}

TEST(EigenvalueErrorTest, NegativeForEnergyBelowExact)
{
    double exact = tise::analyticHydrogenEnergy(1, 0);
    EXPECT_LT(tise::eigenvalueError(exact - 0.001, 1, 0), 0.0);
}

TEST(EigenvalueErrorTest, KnownDifference)
{
    // computed = -0.4, exact for n=1,L=0 = -0.5 => error = 0.1
    EXPECT_NEAR(tise::eigenvalueError(-0.4, 1, 0), 0.1, 1e-15);
}

// ---------------------------------------------------------------------------
// eigenstateCoefficients
// ---------------------------------------------------------------------------

class EigenstateCoefficientsTest : public ::testing::Test
{
protected:
    // nEn=3, nBSplines=5, one eigenvector per column of 3x3 matrix
    // evec = [1,2,3, 4,5,6, 7,8,9] (column-major: col0={1,2,3}, col1={4,5,6}, ...)
    int nEn = 3;
    int nBSplines = 5;
    std::vector<double> evec = {1,2,3, 4,5,6, 7,8,9};
};

TEST_F(EigenstateCoefficientsTest, SizeIsNBSplines)
{
    auto c = tise::eigenstateCoefficients(evec, 1, nEn, nBSplines);
    EXPECT_EQ(static_cast<int>(c.size()), nBSplines);
}

TEST_F(EigenstateCoefficientsTest, FirstAndLastAreZero)
{
    auto c = tise::eigenstateCoefficients(evec, 1, nEn, nBSplines);
    EXPECT_EQ(c.front(), 0.0);
    EXPECT_EQ(c.back(),  0.0);
}

TEST_F(EigenstateCoefficientsTest, InteriorMatchesEvecColumn1)
{
    // Column 1 (1-based iEn=1) is evec[0..2] = {1,2,3}
    auto c = tise::eigenstateCoefficients(evec, 1, nEn, nBSplines);
    EXPECT_EQ(c[1], 1.0);
    EXPECT_EQ(c[2], 2.0);
    EXPECT_EQ(c[3], 3.0);
}

TEST_F(EigenstateCoefficientsTest, InteriorMatchesEvecColumn2)
{
    // Column 2 (1-based iEn=2) is evec[3..5] = {4,5,6}
    auto c = tise::eigenstateCoefficients(evec, 2, nEn, nBSplines);
    EXPECT_EQ(c[1], 4.0);
    EXPECT_EQ(c[2], 5.0);
    EXPECT_EQ(c[3], 6.0);
}

// ---------------------------------------------------------------------------
// fillBandedMatrices
// ---------------------------------------------------------------------------

class FillBandedMatricesTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> grid(nNodes);
        for (int i = 0; i < nNodes; ++i)
            grid[i] = rMin + (rMax - rMin) * i / (nNodes - 1);
        ASSERT_EQ(bs.init(nNodes, order, grid), 0);
        nBSplines = bs.getNBSplines();
        nEn       = nBSplines - 2;

        // Generic two-piece potential (unrelated to hydrogen physics), used only to
        // exercise fillBandedMatrices' multi-domain string/muparser mechanics.
        // Continuous at the split (x=5: "x"->5, "x*x-20.0"->5).
        std::map<std::string, std::string> potential = {
            {"[0.1, 5.0)",  "x"},
            {"[5.0, 10.0]", "x * x - 20.0"}
        };
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L, potential);
        Hmat = H;
        Smat = S;
    }

    bspline::BSpline bs;
    int nNodes = 11, order = 4, L = 0;
    double rMin = 0.1, rMax = 10.0; // avoid x=0 for the potential
    int nBSplines = 0, nEn = 0;
    std::vector<double> Hmat, Smat;

    // Retrieve element (row,col) from a symmetric banded matrix stored in
    // LAPACK 'U' format with leading dimension `ld = order`.
    double bandElem(const std::vector<double> &mat, int row, int col) const
    {
        // row, col: 1-based indices into the nEn x nEn matrix
        if (row > col) std::swap(row, col); // use upper triangle
        if (col - row >= order) return 0.0;
        int brow = row + order - col; // 1..order
        int bcol = col;               // 1..nEn
        return mat[(brow - 1) + (bcol - 1) * order];
    }
};

TEST_F(FillBandedMatricesTest, MatrixSizeIsOrderTimesNEn)
{
    EXPECT_EQ(static_cast<int>(Hmat.size()), order * (nEn + 1));
    EXPECT_EQ(static_cast<int>(Smat.size()), order * (nEn + 1));
}

TEST_F(FillBandedMatricesTest, OverlapDiagonalIsPositive)
{
    for (int i = 1; i <= nEn; ++i)
        EXPECT_GT(bandElem(Smat, i, i), 0.0)
            << "S diagonal not positive at i=" << i;
}

TEST_F(FillBandedMatricesTest, OverlapSymmetric)
{
    for (int i = 1; i <= nEn; ++i)
        for (int j = i; j <= std::min(i + order - 1, nEn); ++j)
            EXPECT_NEAR(bandElem(Smat, i, j), bandElem(Smat, j, i), 1e-11)
                << "S not symmetric at (" << i << "," << j << ")";
}

TEST_F(FillBandedMatricesTest, HamiltonianSymmetric)
{
    for (int i = 1; i <= nEn; ++i)
        for (int j = i; j <= std::min(i + order - 1, nEn); ++j)
            EXPECT_NEAR(bandElem(Hmat, i, j), bandElem(Hmat, j, i), 1e-11)
                << "H not symmetric at (" << i << "," << j << ")";
}

TEST_F(FillBandedMatricesTest, OverlapMatchesDirectIntegral)
{
    auto unity = [](double, const double *) { return 1.0; };
    // S(2,2) in 1-based B-spline is S(1,1) in 1-based nEn matrix (iBs=2)
    double direct = bs.integral(unity, 2, 2);
    EXPECT_NEAR(bandElem(Smat, 1, 1), direct, 1e-11);
}

// ---------------------------------------------------------------------------
// fillBandedMatrices — regression test for the exact radial potential
// ---------------------------------------------------------------------------

class FillBandedMatricesRadialPotentialTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> grid(nNodes);
        for (int i = 0; i < nNodes; ++i)
            grid[i] = rMin + (rMax - rMin) * i / (nNodes - 1);
        ASSERT_EQ(bs.init(nNodes, order, grid), 0);
        nBSplines = bs.getNBSplines();
        nEn       = nBSplines - 2;

        // Muparser encoding of tise::radialPotential(x, L).
        std::map<std::string, std::string> potential = {
            {"(-inf, inf)", std::to_string(L) + " * (" + std::to_string(L) + " + 1.0) / (2.0 * x * x) - 1.0 / x"}
        };
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L, potential);
        Hmat = H;
        Smat = S;
    }

    bspline::BSpline bs;
    int nNodes = 11, order = 4, L = 1; // L=1 so the centrifugal term is nontrivial
    double rMin = 0.1, rMax = 10.0;
    int nBSplines = 0, nEn = 0;
    std::vector<double> Hmat, Smat;

    double bandElem(const std::vector<double> &mat, int row, int col) const
    {
        if (row > col) std::swap(row, col);
        if (col - row >= order) return 0.0;
        int brow = row + order - col;
        int bcol = col;
        return mat[(brow - 1) + (bcol - 1) * order];
    }
};

TEST_F(FillBandedMatricesRadialPotentialTest, HamiltonianMatchesDirectRadialPotentialIntegral)
{
    auto unity = [](double, const double *) { return 1.0; };
    auto directPotential = [this](double x, const double *) {
        return tise::radialPotential(x, L);
    };

    // nEn-index band entry (i,j) <-> B-spline indices (iBs1,iBs2)=(i+1,j+1),
    // per fillBandedMatrices' index algebra (same relation OverlapMatchesDirectIntegral
    // uses implicitly for i=j=1 -> bs.integral(unity,2,2)).
    for (int i = 1; i <= nEn; ++i)
    {
        for (int j = i; j <= std::min(i + order - 1, nEn); ++j)
        {
            int iBs1 = i + 1;
            int iBs2 = j + 1;
            double kinetic   = bs.integral(unity, iBs1, iBs2, 1, 1) / 2.0;
            double potential = bs.integral(directPotential, iBs1, iBs2);
            double expected  = kinetic + potential;
            EXPECT_NEAR(bandElem(Hmat, i, j), expected, 1e-11)
                << "H mismatch at (" << i << "," << j << ")";
        }
    }
}

// ---------------------------------------------------------------------------
// solveGeneralizedEigenproblem
// ---------------------------------------------------------------------------

class SolveEigenTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> grid(nNodes);
        for (int i = 0; i < nNodes; ++i)
            grid[i] = rMin + (rMax - rMin) * i / (nNodes - 1);
        bspline::BSpline bs;
        ASSERT_EQ(bs.init(nNodes, order, grid), 0);
        int nBs = bs.getNBSplines();
        nEn = nBs - 2;

        // Split into two domains with the identical formula, so this test exercises
        // multi-piece domain matching while the analytic-eigenvalue checks below
        // remain valid unchanged.
        double rMid = rMin + (rMax - rMin) / 2.0; // = 30.0, exactly a grid node
        std::string expr = std::to_string(L) + " * (" + std::to_string(L) + " + 1.0) / (2.0 * x * x) - 1.0 / x";
        std::map<std::string, std::string> potential = {
            {"[" + std::to_string(rMin) + ", " + std::to_string(rMid) + ")", expr},
            {"[" + std::to_string(rMid) + ", " + std::to_string(rMax) + "]", expr}
        };
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L, potential);
        result = tise::solveGeneralizedEigenproblem(H, S, nEn, order);
    }

    int nNodes = 31, order = 8, L = 0;
    double rMin = 0.0, rMax = 60.0;
    int nEn = 0;
    tise::EigenResult result;
};

TEST_F(SolveEigenTest, DimensionIsNEn)
{
    EXPECT_EQ(result.dim, nEn);
    EXPECT_EQ(static_cast<int>(result.values.size()),  nEn);
    EXPECT_EQ(static_cast<int>(result.vectors.size()), nEn * nEn);
}

TEST_F(SolveEigenTest, EigenvaluesAscending)
{
    for (int i = 0; i + 1 < result.dim; ++i)
        EXPECT_LE(result.values[i], result.values[i + 1])
            << "Eigenvalues not ascending at i=" << i;
}

TEST_F(SolveEigenTest, GroundStateEnergyL0)
{
    // n=1, L=0 exact = -0.5 a.u.; 31-node/order-8 basis is coarse, use 1e-4
    double exact = tise::analyticHydrogenEnergy(1, 0);
    EXPECT_NEAR(result.values[0], exact, 1e-4);
}

TEST_F(SolveEigenTest, FirstFewEigenvaluesMatchAnalytic)
{
    for (int n = 1; n <= 4; ++n)
    {
        double exact = tise::analyticHydrogenEnergy(n, L);
        EXPECT_NEAR(result.values[n - 1], exact, 1e-4)
            << "Eigenvalue mismatch at n=" << n;
    }
}

// ---------------------------------------------------------------------------
// precomputeBoundaryCoupling
//
// Per docs/planning/tise-task-breakdown.md §3 "B1. Precompute boundary-
// coupling elements": validate <phi_n|B_N> and <phi_n|H|B_N> against a
// directly-computed (brute-force bs.integral call, no shortcuts) reference,
// for a small basis size. B_N here is the true "last B-spline... normally
// dropped to enforce psi(R)=0 for bound states" (1-based bs index nBSplines,
// i.e. nEn+2) -- NOT B_{nEn+1}, which is already the last spline spanning the
// confined nEn-dimensional bound-state basis itself. Using B_{nEn+1} as B_N
// would put B_N inside span{phi_n}, which forces <phi_n|H|B_N> = E_n<phi_n|B_N>
// identically (since H*phi_n = E_n*S*phi_n), collapsing buildContinuumState's
// c_n(E) to a constant independent of E -- see git history for this fix.
// ---------------------------------------------------------------------------

class PrecomputeBoundaryCouplingTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> grid(nNodes);
        for (int i = 0; i < nNodes; ++i)
            grid[i] = rMin + (rMax - rMin) * i / (nNodes - 1);
        ASSERT_EQ(bs.init(nNodes, order, grid), 0);
        nBSplines = bs.getNBSplines();
        nEn       = nBSplines - 2;

        // Muparser encoding of tise::radialPotential(x, L), as in
        // FillBandedMatricesRadialPotentialTest.
        std::map<std::string, std::string> potential = {
            {"(-inf, inf)", std::to_string(L) + " * (" + std::to_string(L) + " + 1.0) / (2.0 * x * x) - 1.0 / x"}
        };
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L, potential);
        Hmat = H;
        Smat = S;
        eigen = tise::solveGeneralizedEigenproblem(H, S, nEn, order);
    }

    bspline::BSpline bs;
    int nNodes = 11, order = 4, L = 1; // small basis; L=1 so H is non-trivial
    double rMin = 0.1, rMax = 10.0;
    int nBSplines = 0, nEn = 0;
    std::vector<double> Hmat, Smat;
    tise::EigenResult eigen;
};

TEST_F(PrecomputeBoundaryCouplingTest, ReturnsOneCoefficientPerEigenstate)
{
    auto [coeffs1, coeffs2] = tise::precomputeBoundaryCoupling(order, nEn, Hmat, Smat, eigen);
    EXPECT_EQ(static_cast<int>(coeffs1.size()), nEn);
    EXPECT_EQ(static_cast<int>(coeffs2.size()), nEn);
}

TEST_F(PrecomputeBoundaryCouplingTest, MatchesBruteForceIntegrals)
{
    auto [coeffs1, coeffs2] = tise::precomputeBoundaryCoupling(order, nEn, Hmat, Smat, eigen);
    ASSERT_EQ(static_cast<int>(coeffs1.size()), nEn);
    ASSERT_EQ(static_cast<int>(coeffs2.size()), nEn);

    auto unity  = [](double, const double *) { return 1.0; };
    auto potFun = [this](double x, const double *) { return tise::radialPotential(x, L); };

    // B_N: the true dropped last B-spline (1-based bs index nBSplines),
    // outside the confined nEn-dimensional bound-state basis.
    int iBsN = nBSplines;

    // Brute-force <B_i|B_N> and <B_i|H|B_N> for each confined B-spline
    // i = 2..nEn+1 (0-based k = i-2), independent of fillBandedMatrices'
    // banded storage / precomputeBoundaryCoupling's HB_N extraction.
    std::vector<double> S_row(nEn), H_row(nEn);
    for (int k = 0; k < nEn; ++k)
    {
        int iBs = k + 2;
        S_row[k] = bs.integral(unity, iBs, iBsN);
        double kinetic   = bs.integral(unity, iBs, iBsN, 1, 1) / 2.0;
        double potential = bs.integral(potFun, iBs, iBsN);
        H_row[k] = kinetic + potential;
    }

    // phi_n = sum_k c_{k,n} |B_{k+2}>, so <phi_n|B_N> = sum_k c_{k,n} <B_{k+2}|B_N>
    // and <phi_n|H|B_N> = sum_k c_{k,n} <B_{k+2}|H|B_N>. eigen.vectors is
    // column-major with leading dimension eigen.ldz.
    for (int n = 0; n < nEn; ++n)
    {
        double expectedOverlap = 0.0, expectedHCoupling = 0.0;
        for (int k = 0; k < nEn; ++k)
        {
            double c = eigen.vectors[n * eigen.ldz + k];
            expectedOverlap   += c * S_row[k];
            expectedHCoupling += c * H_row[k];
        }
        EXPECT_NEAR(coeffs2[n], expectedOverlap,   1e-9) << "<phi_n|B_N> mismatch at n=" << n;
        EXPECT_NEAR(coeffs1[n], expectedHCoupling, 1e-9) << "<phi_n|H|B_N> mismatch at n=" << n;
    }
}

// ---------------------------------------------------------------------------
// buildContinuumState
//
// Per docs/planning/tise-task-breakdown.md §3 "B2. Energy-grid loop and
// closed-form coefficients": for each energy on a grid, compute
// c_n = (<phi_n|H|B_N> - E<phi_n|B_N>) / (E - E_n) and
// |psi_bar_E> = sum_n |phi_n> c_n + |B_N>.
// ---------------------------------------------------------------------------

class BuildContinuumStateTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> grid(nNodes);
        for (int i = 0; i < nNodes; ++i)
            grid[i] = rMin + (rMax - rMin) * i / (nNodes - 1);
        ASSERT_EQ(bs.init(nNodes, order, grid), 0);
        nBSplines = bs.getNBSplines();
        nEn       = nBSplines - 2;

        std::map<std::string, std::string> potential = {
            {"(-inf, inf)", std::to_string(L) + " * (" + std::to_string(L) + " + 1.0) / (2.0 * x * x) - 1.0 / x"}
        };
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L, potential);
        Hmat = H;
        Smat = S;
        eigen = tise::solveGeneralizedEigenproblem(H, S, nEn, order);

        // Energy grid well above the bound-state spectrum (eigen.values are
        // all < 2 for this basis), so E - E_n never vanishes.
        Egrid = {2.0, 5.0, 10.0};
    }

    bspline::BSpline bs;
    int nNodes = 11, order = 4, L = 1;
    double rMin = 0.1, rMax = 10.0;
    int nBSplines = 0, nEn = 0;
    std::vector<double> Hmat, Smat;
    tise::EigenResult eigen;
    std::vector<double> Egrid;
};

TEST_F(BuildContinuumStateTest, ReturnsOneStatePerEnergyGridPoint)
{
    auto states = tise::buildContinuumState(order, nEn, Hmat, Smat, eigen, Egrid);
    EXPECT_EQ(static_cast<int>(states.size()), static_cast<int>(Egrid.size()));
}

TEST_F(BuildContinuumStateTest, EachStateHasNEnPlusOneCoefficients)
{
    auto states = tise::buildContinuumState(order, nEn, Hmat, Smat, eigen, Egrid);
    for (const auto &state : states)
        EXPECT_EQ(static_cast<int>(state.size()), nEn + 1);
}

TEST_F(BuildContinuumStateTest, LastCoefficientIsAlwaysOne)
{
    // The B_N term's coefficient in |psi_bar_E> = sum_n |phi_n> c_n + |B_N>
    // is exactly 1, for every energy on the grid.
    auto states = tise::buildContinuumState(order, nEn, Hmat, Smat, eigen, Egrid);
    for (size_t e = 0; e < states.size(); ++e)
        EXPECT_DOUBLE_EQ(states[e][nEn], 1.0) << "B_N coefficient wrong at E_idx=" << e;
}

TEST_F(BuildContinuumStateTest, CoefficientsVaryWithEnergy)
{
    // Regression guard: c_n(E) must actually depend on E. (It collapses to a
    // constant, independent of E, if B_N is mistakenly drawn from inside
    // span{phi_n} rather than the true dropped last B-spline -- see
    // PrecomputeBoundaryCouplingTest above.)
    auto states = tise::buildContinuumState(order, nEn, Hmat, Smat, eigen, Egrid);
    ASSERT_GE(states.size(), 2u);

    bool anyDifferent = false;
    for (int n = 0; n < nEn; ++n)
        if (std::abs(states[0][n] - states[1][n]) > 1e-9)
            anyDifferent = true;
    EXPECT_TRUE(anyDifferent) << "coefficients identical across different energies";
}

TEST_F(BuildContinuumStateTest, SatisfiesDefiningEigenrelation)
{
    // The doc's B2 "Done when" criterion: <phi_n|(E-H)|psi_bar_E> ~ 0 for
    // every confined n, at each sampled energy. |psi_bar_E> = sum_n c_n|phi_n>
    // + |B_N>, and {phi_n} is S-orthonormal with H phi_n = E_n S phi_n (the
    // defining property of solveGeneralizedEigenproblem's output), so this
    // reduces to E*(c_n + <phi_n|B_N>) - (c_n*E_n + <phi_n|H|B_N>) ~ 0.
    // <phi_n|B_N> and <phi_n|H|B_N> are precomputeBoundaryCoupling's
    // coeffs2/coeffs1, independently validated against brute-force
    // bs.integral() calls by PrecomputeBoundaryCouplingTest above -- so this
    // is not simply re-checking the algebra that defined c_n.
    auto states = tise::buildContinuumState(order, nEn, Hmat, Smat, eigen, Egrid);
    auto [coeffs1, coeffs2] = tise::precomputeBoundaryCoupling(order, nEn, Hmat, Smat, eigen);

    for (size_t eIdx = 0; eIdx < Egrid.size(); ++eIdx)
    {
        double E = Egrid[eIdx];
        for (int n = 0; n < nEn; ++n)
        {
            double c = states[eIdx][n];
            double lhs = E * (c + coeffs2[n]) - (c * eigen.values[n] + coeffs1[n]);
            EXPECT_NEAR(lhs, 0.0, 1e-9)
                << "<phi_n|(E-H)|psi_bar_E> nonzero at n=" << n << ", E=" << E;
        }
    }
}

// ---------------------------------------------------------------------------
// writeEigenstate
// ---------------------------------------------------------------------------

class WriteEigenstateTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> grid(nNodes);
        for (int i = 0; i < nNodes; ++i)
            grid[i] = i * 1.0 / (nNodes - 1);
        ASSERT_EQ(bs.init(nNodes, order, grid), 0);
        nBSplines = bs.getNBSplines();
        coeffs.assign(nBSplines, 0.0);
        coeffs[nBSplines / 2] = 1.0; // arbitrary non-trivial wavefunction
    }

    bspline::BSpline bs;
    int nNodes = 11, order = 4, nBSplines = 0;
    std::vector<double> coeffs;
};

TEST_F(WriteEigenstateTest, LineCountEqualsNpts)
{
    std::ostringstream ss;
    tise::writeEigenstate(ss, bs, coeffs, 50, 0.0, 1.0);
    std::string out = ss.str();
    int lines = static_cast<int>(std::count(out.begin(), out.end(), '\n'));
    EXPECT_EQ(lines, 50);
}

TEST_F(WriteEigenstateTest, XValuesAreMonotone)
{
    std::ostringstream ss;
    tise::writeEigenstate(ss, bs, coeffs, 10, 0.0, 1.0);
    std::istringstream in(ss.str());
    double prevX = -1.0, x, y;
    while (in >> x >> y)
    {
        EXPECT_GT(x, prevX);
        prevX = x;
    }
}

TEST_F(WriteEigenstateTest, FirstXIsRMin)
{
    std::ostringstream ss;
    tise::writeEigenstate(ss, bs, coeffs, 5, 0.1, 0.9);
    std::istringstream in(ss.str());
    double x, y;
    in >> x >> y;
    EXPECT_NEAR(x, 0.1, 1e-14);
}

TEST_F(WriteEigenstateTest, LastXIsRMax)
{
    const int npts = 5;
    std::ostringstream ss;
    tise::writeEigenstate(ss, bs, coeffs, npts, 0.1, 0.9);
    std::istringstream in(ss.str());
    double x = 0, y = 0, lastX = 0;
    while (in >> x >> y) lastX = x;
    EXPECT_NEAR(lastX, 0.9, 1e-14);
}

// ---------------------------------------------------------------------------
// matchAsymptotic — phase shift vs. the exact spherical square-well solution
//
// For an L=0 (s-wave) attractive square well V(r) = -V0 for r < a, 0 for
// r >= a, the continuum phase shift at energy E > 0 has a closed form (finite
// spherical well scattering):
//   k     = sqrt(E)          exterior wavenumber
//   kappa = sqrt(E + V0)     interior wavenumber
//   delta = -k*a + atan[ (k/kappa) * tan(kappa*a) ] + n*pi, for integer n.
// delta is only physically defined mod pi (sin(kr+delta) and
// sin(kr+delta+pi) are the same scattering state up to an overall sign), so
// comparisons below reduce both sides into (-pi/2, pi/2] before comparing.
// ---------------------------------------------------------------------------

static double wrapPhaseModPi(double delta)
{
    return delta - M_PI * std::round(delta / M_PI);
}

static double squareWellPhaseShift(double E, double V0, double a)
{
    double k = std::sqrt(2*E);
    double kappa = std::sqrt(2*(E + V0));
    return -k * a + std::atan((k / kappa) * std::tan(kappa * a));
}

class SquareWellPhaseShiftTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> gridPts(nNodes);
        for (int i = 0; i < nNodes; ++i)
            gridPts[i] = rMin + (rMax - rMin) * i / (nNodes - 1);
        ASSERT_EQ(bs.init(nNodes, order, gridPts), 0);
        nBSplines = bs.getNBSplines();
        nEn       = nBSplines - 2;

        // Attractive square well of depth V0 and radius a; a coincides with a
        // grid node (spacing (rMax-rMin)/(nNodes-1) = 0.5) so a B-spline knot
        // sits exactly at the potential's discontinuity instead of smoothing
        // over it.
        std::map<std::string, std::string> potential = {
            {"[0, " + std::to_string(a) + ")", "-" + std::to_string(V0)},
            {"[" + std::to_string(a) + ", " + std::to_string(rMax) + "]", "0.0"}
        };
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L, potential);
        eigen  = tise::solveGeneralizedEigenproblem(H, S, nEn, order);
        states = tise::buildContinuumState(order, nEn, H, S, eigen, energyGrid);
    }

    bspline::BSpline bs;
    int nNodes = 41, order = 6, L = 0;
    double rMin = 0.0, rMax = 20.0;
    double V0 = 2.0, a = 0.5;
    int nBSplines = 0, nEn = 0;
    tise::EigenResult eigen;
    std::vector<double> energyGrid = {0.5};
    std::vector<std::vector<double>> states;
};

TEST_F(SquareWellPhaseShiftTest, MatchesAnalyticSquareWellFormula)
{
    double R = 15.0; // well outside the well (a=0.5), well inside rMax=20
    auto ar = tise::matchAsymptotic(bs, states, eigen, energyGrid, R);

    double expected = squareWellPhaseShift(energyGrid[0], V0, a);
    EXPECT_NEAR(wrapPhaseModPi(ar.delta[0]), wrapPhaseModPi(expected), 5e-3)
        << "computed delta=" << ar.delta[0] << " expected=" << expected;
}

TEST_F(SquareWellPhaseShiftTest, PhaseShiftIndependentOfMatchingRadius)
{
    // Physically, once R is outside the well, the continuum wavefunction is
    // an exact standing wave A*sin(kr+delta); delta extracted at different R
    // (mod pi) must therefore agree with itself, independent of the
    // analytic-formula tolerance used above.
    std::vector<double> Rs = {10.0, 12.0, 15.0, 18.0};
    double deltaAtR0 = wrapPhaseModPi(
        tise::matchAsymptotic(bs, states, eigen, energyGrid, Rs[0]).delta[0]);

    for (double R : Rs)
    {
        double delta = wrapPhaseModPi(
            tise::matchAsymptotic(bs, states, eigen, energyGrid, R).delta[0]);
        EXPECT_NEAR(delta, deltaAtR0, 1e-3)
            << "phase shift not independent of matching radius R=" << R;
    }
}
