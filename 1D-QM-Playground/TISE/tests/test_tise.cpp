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
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L);
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
    EXPECT_EQ(static_cast<int>(Hmat.size()), order * nEn);
    EXPECT_EQ(static_cast<int>(Smat.size()), order * nEn);
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
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L);
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
