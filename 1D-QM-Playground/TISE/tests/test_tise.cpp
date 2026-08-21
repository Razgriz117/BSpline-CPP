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
// inInterval
// ---------------------------------------------------------------------------

TEST(InIntervalTest, ClosedBothEnds)
{
    EXPECT_TRUE(tise::inInterval(0.0, "[0, 5]"));
    EXPECT_TRUE(tise::inInterval(5.0, "[0, 5]"));
    EXPECT_TRUE(tise::inInterval(2.5, "[0, 5]"));
    EXPECT_FALSE(tise::inInterval(-0.1, "[0, 5]"));
    EXPECT_FALSE(tise::inInterval(5.1, "[0, 5]"));
}

TEST(InIntervalTest, HalfOpenVariants)
{
    EXPECT_FALSE(tise::inInterval(0.0, "(0, 5)"));
    EXPECT_FALSE(tise::inInterval(5.0, "(0, 5)"));
    EXPECT_TRUE(tise::inInterval(2.5, "(0, 5)"));

    EXPECT_TRUE(tise::inInterval(0.0, "[0, 5)"));
    EXPECT_FALSE(tise::inInterval(5.0, "[0, 5)"));

    EXPECT_FALSE(tise::inInterval(0.0, "(0, 5]"));
    EXPECT_TRUE(tise::inInterval(5.0, "(0, 5]"));
}

TEST(InIntervalTest, InfinityBounds)
{
    EXPECT_TRUE(tise::inInterval(1e10, "[0, inf)"));
    EXPECT_FALSE(tise::inInterval(-1.0, "[0, inf)"));
    EXPECT_TRUE(tise::inInterval(-1e10, "(-inf, 0]"));
    EXPECT_FALSE(tise::inInterval(0.1, "(-inf, 0]"));
    // Case-insensitive "infinity" spelling.
    EXPECT_TRUE(tise::inInterval(1e10, "[0, infinity)"));
}

TEST(InIntervalTest, ThrowsOnMalformedString)
{
    EXPECT_THROW(tise::inInterval(1.0, "not an interval"), std::runtime_error);
}

// ---------------------------------------------------------------------------
// evaluateFunction
// ---------------------------------------------------------------------------

TEST(EvaluateFunctionTest, PicksCorrectPiece)
{
    std::map<std::string, std::string> potential = {
        {"[0, 5)", "x * x"},
        {"[5, 10]", "10 - x"}
    };
    EXPECT_NEAR(tise::evaluateFunction(potential, 3.0), 9.0, 1e-12);
    EXPECT_NEAR(tise::evaluateFunction(potential, 7.0), 3.0, 1e-12);
}

TEST(EvaluateFunctionTest, ThrowsWhenXUncovered)
{
    std::map<std::string, std::string> potential = {
        {"[0, 5)", "x"}
    };
    EXPECT_THROW(tise::evaluateFunction(potential, 10.0), std::runtime_error);
}

// ---------------------------------------------------------------------------
// classifySequenceConvergence
// ---------------------------------------------------------------------------

TEST(ClassifySequenceConvergenceTest, DetectsFlatSequence)
{
    std::vector<double> V(16, 10.0); // exactly constant
    auto fit = tise::classifySequenceConvergence(V, 4.0);
    EXPECT_FALSE(fit.isDivergent);
    EXPECT_TRUE(fit.isFlat);
    EXPECT_NEAR(fit.fittedLimit, 10.0, 1e-12);
}

TEST(ClassifySequenceConvergenceTest, DetectsDivergentSequence)
{
    std::vector<double> V(16);
    for (int k = 0; k < 16; ++k)
        V[k] = std::pow(16.0, k); // differences grow geometrically, never shrink
    auto fit = tise::classifySequenceConvergence(V, 4.0);
    EXPECT_TRUE(fit.isDivergent);
    EXPECT_FALSE(fit.isFlat);
}

TEST(ClassifySequenceConvergenceTest, FitsCoulombPowerLaw)
{
    // V[k] = 3 * ratio^(-k): a pure p=1 power-law decay to a zero limit.
    std::vector<double> V(16);
    for (int k = 0; k < 16; ++k)
        V[k] = 3.0 * std::pow(4.0, -k);
    auto fit = tise::classifySequenceConvergence(V, 4.0);
    EXPECT_FALSE(fit.isDivergent);
    EXPECT_FALSE(fit.isFlat);
    EXPECT_NEAR(fit.powerLawExponent, 1.0, 1e-6);
    EXPECT_NEAR(fit.fittedLimit, 0.0, 1e-9);
}

TEST(ClassifySequenceConvergenceTest, FitsArbitraryPowerLaw)
{
    // V[k] = 3 * ratio^(-1.5k): a pure p=1.5 power-law decay to a zero limit.
    std::vector<double> V(16);
    for (int k = 0; k < 16; ++k)
        V[k] = 3.0 * std::pow(4.0, -1.5 * k);
    auto fit = tise::classifySequenceConvergence(V, 4.0);
    EXPECT_FALSE(fit.isDivergent);
    EXPECT_FALSE(fit.isFlat);
    EXPECT_NEAR(fit.powerLawExponent, 1.5, 1e-6);
    EXPECT_NEAR(fit.fittedLimit, 0.0, 1e-9);
}

// ---------------------------------------------------------------------------
// case3WindowFunction / evaluateWindowedPotential
// (docs/planning/boundary-condition-case-3-smoothing.md)
// ---------------------------------------------------------------------------

TEST(Case3WindowFunctionTest, EqualsOneWellInsideBoundary)
{
    EXPECT_NEAR(tise::case3WindowFunction(5.0, 10.0, 2.0, tise::DomainSide::Right), 1.0, 1e-14);
    EXPECT_NEAR(tise::case3WindowFunction(8.0, 10.0, 2.0, tise::DomainSide::Right), 1.0, 1e-14);
}

TEST(Case3WindowFunctionTest, EqualsZeroAtAndBeyondBoundary)
{
    EXPECT_NEAR(tise::case3WindowFunction(10.0, 10.0, 2.0, tise::DomainSide::Right), 0.0, 1e-14);
    EXPECT_NEAR(tise::case3WindowFunction(12.0, 10.0, 2.0, tise::DomainSide::Right), 0.0, 1e-14);
}

TEST(Case3WindowFunctionTest, IsContinuousAtBothTransitionEdges)
{
    const double R = 10.0, delta = 2.0;
    const double eps = 1e-6;
    EXPECT_NEAR(tise::case3WindowFunction(R - delta - eps, R, delta, tise::DomainSide::Right),
                tise::case3WindowFunction(R - delta + eps, R, delta, tise::DomainSide::Right), 1e-4);
    EXPECT_NEAR(tise::case3WindowFunction(R - eps, R, delta, tise::DomainSide::Right),
                tise::case3WindowFunction(R + eps, R, delta, tise::DomainSide::Right), 1e-4);
}

TEST(Case3WindowFunctionTest, DerivativeVanishesAtBothTransitionEdges)
{
    const double R = 10.0, delta = 2.0;
    const double h = 1e-5;
    auto centralDiff = [&](double x) {
        return (tise::case3WindowFunction(x + h, R, delta, tise::DomainSide::Right) -
                tise::case3WindowFunction(x - h, R, delta, tise::DomainSide::Right)) / (2.0 * h);
    };
    EXPECT_NEAR(centralDiff(R - delta), 0.0, 1e-3);
    EXPECT_NEAR(centralDiff(R), 0.0, 1e-3);
}

TEST(Case3WindowFunctionTest, MonotonicWithinTransition)
{
    const double R = 10.0, delta = 2.0;
    double prev = tise::case3WindowFunction(R - delta, R, delta, tise::DomainSide::Right);
    for (int i = 1; i <= 10; ++i)
    {
        double x = (R - delta) + delta * i / 10.0;
        double w = tise::case3WindowFunction(x, R, delta, tise::DomainSide::Right);
        EXPECT_LE(w, prev + 1e-12);
        prev = w;
    }
}

TEST(Case3WindowFunctionTest, LeftSideMirrorsRightSide)
{
    const double R = 10.0, delta = 2.0;
    // Left side: interior (trusted) region is x >= R+delta; wall is at x <= R.
    EXPECT_NEAR(tise::case3WindowFunction(R + delta, R, delta, tise::DomainSide::Left), 1.0, 1e-14);
    EXPECT_NEAR(tise::case3WindowFunction(R + 5.0, R, delta, tise::DomainSide::Left), 1.0, 1e-14);
    EXPECT_NEAR(tise::case3WindowFunction(R, R, delta, tise::DomainSide::Left), 0.0, 1e-14);
    EXPECT_NEAR(tise::case3WindowFunction(R - 5.0, R, delta, tise::DomainSide::Left), 0.0, 1e-14);
}

TEST(EvaluateWindowedPotentialTest, MatchesRawPotentialWellInsideBoundary)
{
    std::map<std::string, std::string> potential = {{"(-inf, inf)", "x*x"}};
    double result = tise::evaluateWindowedPotential(potential, 5.0, 10.0, 2.0, tise::DomainSide::Right);
    EXPECT_NEAR(result, tise::evaluateFunction(potential, 5.0), 1e-12);
}

TEST(EvaluateWindowedPotentialTest, VanishesAtAndBeyondBoundary)
{
    std::map<std::string, std::string> potential = {{"(-inf, inf)", "x*x"}};
    EXPECT_NEAR(tise::evaluateWindowedPotential(potential, 10.0, 10.0, 2.0, tise::DomainSide::Right), 0.0, 1e-12);
    EXPECT_NEAR(tise::evaluateWindowedPotential(potential, 12.0, 10.0, 2.0, tise::DomainSide::Right), 0.0, 1e-12);
}

TEST(EvaluateWindowedPotentialTest, MatchesDirectWindowFunctionTimesRawPotential)
{
    std::map<std::string, std::string> potential = {{"(-inf, inf)", "x*x"}};
    double x = 9.0, R = 10.0, delta = 2.0;
    double expected = tise::case3WindowFunction(x, R, delta, tise::DomainSide::Right) *
                       tise::evaluateFunction(potential, x);
    EXPECT_NEAR(tise::evaluateWindowedPotential(potential, x, R, delta, tise::DomainSide::Right),
                expected, 1e-12);
}

// ---------------------------------------------------------------------------
// classifyAsymptote
// ---------------------------------------------------------------------------

TEST(ClassifyAsymptoteTest, Case1HardWallForQuadraticGrowth)
{
    std::map<std::string, std::string> potential = {{"[0, inf)", "x*x"}};
    tise::SpatialDomain domain{0.0, 20.0};
    std::ostringstream warn;
    auto result = tise::classifyAsymptote(potential, domain, tise::DomainSide::Right, warn);
    EXPECT_EQ(result.asymptoteCase, tise::AsymptoteCase::HardWall);
}

TEST(ClassifyAsymptoteTest, Case1HardWallForLinearGrowth)
{
    std::map<std::string, std::string> potential = {{"[0, inf)", "x"}};
    tise::SpatialDomain domain{0.0, 20.0};
    std::ostringstream warn;
    auto result = tise::classifyAsymptote(potential, domain, tise::DomainSide::Right, warn);
    EXPECT_EQ(result.asymptoteCase, tise::AsymptoteCase::HardWall);
}

TEST(ClassifyAsymptoteTest, Case2FlatForStepPotential)
{
    std::map<std::string, std::string> potential = {{"[0,5)", "0"}, {"[5, inf)", "10"}};
    tise::SpatialDomain domain{0.0, 20.0};
    std::ostringstream warn;
    auto result = tise::classifyAsymptote(potential, domain, tise::DomainSide::Right, warn);
    EXPECT_EQ(result.asymptoteCase, tise::AsymptoteCase::AnalyticAsymptote);
    EXPECT_EQ(result.subType, tise::AsymptoteSubType::Flat);
    EXPECT_NEAR(result.fittedAsymptoticValue, 10.0, 1e-6);
}

TEST(ClassifyAsymptoteTest, Case2CoulombForInverseR)
{
    std::map<std::string, std::string> potential = {{"(0, inf)", "-1/x"}};
    tise::SpatialDomain domain{0.1, 50.0};
    std::ostringstream warn;
    auto result = tise::classifyAsymptote(potential, domain, tise::DomainSide::Right, warn);
    EXPECT_EQ(result.asymptoteCase, tise::AsymptoteCase::AnalyticAsymptote);
    EXPECT_EQ(result.subType, tise::AsymptoteSubType::Coulomb);
    EXPECT_NEAR(result.powerLawExponent, 1.0, 1e-2);
}

TEST(ClassifyAsymptoteTest, Case3IrregularForPowerLawOneAndHalf)
{
    std::map<std::string, std::string> potential = {{"(0, inf)", "1/x^1.5"}};
    tise::SpatialDomain domain{0.1, 50.0};
    std::ostringstream warn;
    auto result = tise::classifyAsymptote(potential, domain, tise::DomainSide::Right, warn);
    EXPECT_EQ(result.asymptoteCase, tise::AsymptoteCase::Irregular);
    EXPECT_NEAR(result.powerLawExponent, 1.5, 1e-2);
    EXPECT_TRUE(result.warningEmitted);
    EXPECT_GT(result.recommendedTransitionWidth, 0.0);
}

TEST(ClassifyAsymptoteTest, Case3WarningTextMentionsTaperingNotDiscontinuity)
{
    std::map<std::string, std::string> potential = {{"(0, inf)", "1/x^1.5"}};
    tise::SpatialDomain domain{0.1, 50.0};
    std::ostringstream warn;
    tise::classifyAsymptote(potential, domain, tise::DomainSide::Right, warn);
    std::string text = warn.str();
    EXPECT_NE(text.find("taper"), std::string::npos);
    EXPECT_NE(text.find("approximate"), std::string::npos);
    EXPECT_EQ(text.find("discontinuity"), std::string::npos);
}

TEST(ClassifyAsymptoteTest, NoWarningForCase1AndCase2)
{
    tise::SpatialDomain domain{0.0, 20.0};
    std::ostringstream warn1;
    std::map<std::string, std::string> hardWallPotential = {{"[0, inf)", "x*x"}};
    tise::classifyAsymptote(hardWallPotential, domain, tise::DomainSide::Right, warn1);
    EXPECT_TRUE(warn1.str().empty());

    std::ostringstream warn2;
    std::map<std::string, std::string> flatPotential = {{"[0,5)", "0"}, {"[5, inf)", "10"}};
    tise::classifyAsymptote(flatPotential, domain, tise::DomainSide::Right, warn2);
    EXPECT_TRUE(warn2.str().empty());
}

TEST(ClassifyAsymptoteTest, LeftSideSymmetric)
{
    std::map<std::string, std::string> potential = {{"(-inf, inf)", "x*x"}};
    tise::SpatialDomain domain{-20.0, 0.0};
    std::ostringstream warn;
    auto result = tise::classifyAsymptote(potential, domain, tise::DomainSide::Left, warn);
    EXPECT_EQ(result.asymptoteCase, tise::AsymptoteCase::HardWall);
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
// classifyBoundStates
// ---------------------------------------------------------------------------

TEST(ClassifyBoundStatesTest, AllBoundWhenAllBelowThreshold)
{
    tise::EigenResult r;
    r.values = {-5.0, -3.0, -1.0};
    r.dim = 3;
    auto c = tise::classifyBoundStates(r, 0.0);
    EXPECT_EQ(c.nBound, 3);
    EXPECT_TRUE(c.isBound[0]);
    EXPECT_TRUE(c.isBound[1]);
    EXPECT_TRUE(c.isBound[2]);
}

TEST(ClassifyBoundStatesTest, AllAboveThresholdWhenNoneBelow)
{
    tise::EigenResult r;
    r.values = {0.5, 1.0, 2.0};
    r.dim = 3;
    auto c = tise::classifyBoundStates(r, 0.0);
    EXPECT_EQ(c.nBound, 0);
    for (bool b : c.isBound)
        EXPECT_FALSE(b);
}

TEST(ClassifyBoundStatesTest, MixedAboveAndBelowThreshold)
{
    tise::EigenResult r;
    r.values = {-2.0, -1.0, 0.5, 1.5};
    r.dim = 4;
    auto c = tise::classifyBoundStates(r, 0.0);
    EXPECT_EQ(c.nBound, 2);
    EXPECT_TRUE(c.isBound[0]);
    EXPECT_TRUE(c.isBound[1]);
    EXPECT_FALSE(c.isBound[2]);
    EXPECT_FALSE(c.isBound[3]);
}

TEST(ClassifyBoundStatesTest, ValueExactlyAtThresholdIsNotBound)
{
    tise::EigenResult r;
    r.values = {-1.0, 0.0, 1.0};
    r.dim = 3;
    auto c = tise::classifyBoundStates(r, 0.0);
    EXPECT_TRUE(c.isBound[0]);
    EXPECT_FALSE(c.isBound[1]); // marginal case: strict '<', not bound
    EXPECT_FALSE(c.isBound[2]);
    EXPECT_EQ(c.nBound, 1);
}

class ClassifyBoundStatesHydrogenTest : public ::testing::Test
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

        double rMid = rMin + (rMax - rMin) / 2.0;
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

TEST_F(ClassifyBoundStatesHydrogenTest, GroundStateBoundNBoundLessThanDimAndPrefixStructureHolds)
{
    auto c = tise::classifyBoundStates(result, 0.0);
    ASSERT_FALSE(c.isBound.empty());
    EXPECT_TRUE(c.isBound[0]) << "Ground state (E=-0.5 a.u.) must be bound at threshold 0.0";
    EXPECT_LT(c.nBound, result.dim);

    bool seenFalse = false;
    for (bool b : c.isBound)
    {
        if (!b)
            seenFalse = true;
        else
            EXPECT_FALSE(seenFalse) << "true entry found after a false entry -- not a prefix";
    }
}

class ClassifyBoundStatesSquareWellTest : public ::testing::Test
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

        std::map<std::string, std::string> potential = {
            {"[0, 5)", "-2.0"},
            {"[5, 30]", "0.0"}
        };
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L, potential);
        result = tise::solveGeneralizedEigenproblem(H, S, nEn, order);
    }

    // nNodes=121 (not 31, like the hydrogenic fixture) because this potential's
    // step discontinuity at x=5 isn't yet resolved with degenerate knots (that's
    // task A4's strategic node placement); a uniform grid converges to the
    // analytic odd-parity energies, just more slowly. Empirically verified:
    // nNodes=31 gives ~1e-2 energy error, nNodes=121 gives ~2e-4. See
    // docs/planning/engineer-a-plan-A2.md.
    int nNodes = 121, order = 8, L = 0; // L unused for potential selection; map fully specifies V(x)
    double rMin = 0.0, rMax = 30.0;
    int nEn = 0;
    tise::EigenResult result;
};

TEST_F(ClassifyBoundStatesSquareWellTest, ExactlyThreeBoundStatesBelowZero)
{
    auto c = tise::classifyBoundStates(result, 0.0);
    EXPECT_EQ(c.nBound, 3);
    ASSERT_GT(static_cast<int>(c.isBound.size()), 3);
    EXPECT_TRUE(c.isBound[0]);
    EXPECT_TRUE(c.isBound[1]);
    EXPECT_TRUE(c.isBound[2]);
    EXPECT_FALSE(c.isBound[3]);
}

TEST_F(ClassifyBoundStatesSquareWellTest, BoundEnergiesMatchAnalyticOddParityRoots)
{
    // Analytic odd-parity roots for V0=2.0, w=5.0 (transcendental matching
    // equation k*cot(k*w)=-kappa, cross-checked by independent
    // finite-difference diagonalization -- see docs/planning/engineer-a-plan-A2.md).
    ASSERT_GE(static_cast<int>(result.values.size()), 3);
    EXPECT_NEAR(result.values[0], -1.83728, 1e-3);
    EXPECT_NEAR(result.values[1], -1.35493, 1e-3);
    EXPECT_NEAR(result.values[2], -0.58099, 1e-3);
}

// ---------------------------------------------------------------------------
// checkWellContainment
// ---------------------------------------------------------------------------

class WellContainmentSmallBoxTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> grid(nNodes);
        for (int i = 0; i < nNodes; ++i)
            grid[i] = rMin + (rMax - rMin) * i / (nNodes - 1);
        ASSERT_EQ(bs.init(nNodes, order, grid), 0);
        int nBs = bs.getNBSplines();
        int nEn = nBs - 2;

        double rMid = rMin + (rMax - rMin) / 2.0;
        std::string expr = std::to_string(L) + " * (" + std::to_string(L) + " + 1.0) / (2.0 * x * x) - 1.0 / x";
        std::map<std::string, std::string> potential = {
            {"[" + std::to_string(rMin) + ", " + std::to_string(rMid) + ")", expr},
            {"[" + std::to_string(rMid) + ", " + std::to_string(rMax) + "]", expr}
        };
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L, potential);
        auto result = tise::solveGeneralizedEigenproblem(H, S, nEn, order);
        coeffs = tise::eigenstateCoefficients(result.vectors, 1, nEn, nBs);
    }

    // Deliberately too small: rMax=2.0 is only ~2 decay lengths beyond the
    // hydrogenic ground state's kappa=1, so the ground state "collides" with
    // the wall. |psi'(2.0)| ~ 0.27, ~270x the default tol (see plan doc).
    int nNodes = 31, order = 8, L = 0;
    double rMin = 0.0, rMax = 2.0;
    bspline::BSpline bs;
    std::vector<double> coeffs;
};

TEST_F(WellContainmentSmallBoxTest, FlagsGroundState)
{
    auto c = tise::checkWellContainment(bs, coeffs, rMax);
    EXPECT_TRUE(c.notWellContained);
}

class WellContainmentLargeBoxTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> grid(nNodes);
        for (int i = 0; i < nNodes; ++i)
            grid[i] = rMin + (rMax - rMin) * i / (nNodes - 1);
        ASSERT_EQ(bs.init(nNodes, order, grid), 0);
        int nBs = bs.getNBSplines();
        int nEn = nBs - 2;

        double rMid = rMin + (rMax - rMin) / 2.0;
        std::string expr = std::to_string(L) + " * (" + std::to_string(L) + " + 1.0) / (2.0 * x * x) - 1.0 / x";
        std::map<std::string, std::string> potential = {
            {"[" + std::to_string(rMin) + ", " + std::to_string(rMid) + ")", expr},
            {"[" + std::to_string(rMid) + ", " + std::to_string(rMax) + "]", expr}
        };
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L, potential);
        auto result = tise::solveGeneralizedEigenproblem(H, S, nEn, order);
        coeffs = tise::eigenstateCoefficients(result.vectors, 1, nEn, nBs);
    }

    // Matches SolveEigenTest's converged box -- ~60 decay lengths of margin
    // for the ground state (|psi'(60.0)| ~ 1e-24).
    int nNodes = 31, order = 8, L = 0;
    double rMin = 0.0, rMax = 60.0;
    bspline::BSpline bs;
    std::vector<double> coeffs;
};

TEST_F(WellContainmentLargeBoxTest, DoesNotFlagGroundState)
{
    auto c = tise::checkWellContainment(bs, coeffs, rMax);
    EXPECT_FALSE(c.notWellContained);
}

TEST_F(WellContainmentLargeBoxTest, DerivativeMatchesDirectBsEvalCall)
{
    auto c = tise::checkWellContainment(bs, coeffs, rMax);
    const int n = static_cast<int>(coeffs.size());
    double direct = bs.eval(rMax, coeffs.data(), n, 1);
    EXPECT_NEAR(c.psiPrimeAtBoundary, direct, 1e-15);
}

// Cheap, non-physics fixture (mirrors WriteEigenstateTest's setup) that pins
// checkWellContainment's tol-comparison semantics exactly, independent of
// LAPACK/muparser. Coefficients ramp from 1 (interior) to a hard 0 at both
// ends (mirrors eigenstateCoefficients' zero-padding), robust to exactly
// which B-spline's support reaches the boundary -- WriteEigenstateTest's own
// single-mid-coefficient pattern was checked and does NOT reach x=grid.back()
// for nNodes=11/order=4, so that exact pattern isn't reused here.
class WellContainmentSyntheticTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> grid(nNodes);
        for (int i = 0; i < nNodes; ++i)
            grid[i] = i * 1.0 / (nNodes - 1);
        ASSERT_EQ(bs.init(nNodes, order, grid), 0);
        int nBSplines = bs.getNBSplines();
        coeffs.assign(nBSplines, 1.0);
        coeffs.front() = 0.0;
        coeffs.back() = 0.0;
        xBoundary = grid.back();
    }

    int nNodes = 11, order = 4;
    bspline::BSpline bs;
    std::vector<double> coeffs;
    double xBoundary = 0.0;
};

TEST_F(WellContainmentSyntheticTest, DerivativeExactlyAtToleranceIsNotFlagged)
{
    const int n = static_cast<int>(coeffs.size());
    double raw = bs.eval(xBoundary, coeffs.data(), n, 1);
    ASSERT_NE(raw, 0.0) << "test setup assumption: psi'(xBoundary) must be nonzero";

    auto c = tise::checkWellContainment(bs, coeffs, xBoundary, std::abs(raw));
    EXPECT_NEAR(c.psiPrimeAtBoundary, raw, 1e-15);
    EXPECT_FALSE(c.notWellContained)
        << "exactly at tol is the marginal case: strict '>' required to flag";
}

TEST_F(WellContainmentSyntheticTest, DerivativeJustAboveToleranceIsFlagged)
{
    const int n = static_cast<int>(coeffs.size());
    double raw = bs.eval(xBoundary, coeffs.data(), n, 1);
    ASSERT_NE(raw, 0.0) << "test setup assumption: psi'(xBoundary) must be nonzero";

    auto c = tise::checkWellContainment(bs, coeffs, xBoundary, std::abs(raw) * 0.999999);
    EXPECT_TRUE(c.notWellContained);
}

// ---------------------------------------------------------------------------
// detectPotentialStructure
// ---------------------------------------------------------------------------

TEST(DetectPotentialStructureTest, StepAtBothBarrierEdges)
{
    std::map<std::string, std::string> potential = {
        {"[0,5)", "0"}, {"[5,6]", "10"}, {"(6,10]", "0"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 2u);
    std::sort(joins.begin(), joins.end(),
              [](const tise::DetectedJoin &a, const tise::DetectedJoin &b) { return a.x < b.x; });
    EXPECT_DOUBLE_EQ(joins[0].x, 5.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::Step);
    EXPECT_DOUBLE_EQ(joins[1].x, 6.0);
    EXPECT_EQ(joins[1].type, tise::JoinType::Step);
}

TEST(DetectPotentialStructureTest, StitchedKinkAtCornerJoin)
{
    std::map<std::string, std::string> potential = {
        {"[0,5)", "x"}, {"[5,10]", "10-x"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    EXPECT_DOUBLE_EQ(joins[0].x, 5.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::StitchedKink);
}

TEST(DetectPotentialStructureTest, SingularAtCoulombOrigin)
{
    std::map<std::string, std::string> potential = {
        {"(0, 100]", "-1/x + 1/x^2"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    EXPECT_DOUBLE_EQ(joins[0].x, 0.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::Singular);
}

TEST(DetectPotentialStructureTest, ContinuousJoinNotFlagged)
{
    std::map<std::string, std::string> potential = {
        {"[0,5)", "x*x"}, {"[5,10]", "x*x"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    EXPECT_DOUBLE_EQ(joins[0].x, 5.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::Continuous);
}

TEST(DetectPotentialStructureTest, SingularAtInteriorJoin)
{
    // 1/(1-x) diverges approaching x=1 from the left; exercises
    // classifyJoin's own Singular branch (distinct from the domain-edge
    // singularity path SingularAtCoulombOrigin covers).
    std::map<std::string, std::string> potential = {
        {"[0,1)", "1/(1-x)"}, {"[1,2]", "0"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    EXPECT_DOUBLE_EQ(joins[0].x, 1.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::Singular);
}

TEST(DetectPotentialStructureTest, SingularAtRightDomainEdge)
{
    // 1/(100-x) diverges approaching x=100 from the left; exercises the
    // last-piece (right) global-edge singularity push, distinct from
    // SingularAtCoulombOrigin's left-edge case.
    std::map<std::string, std::string> potential = {
        {"[0, 100)", "1/(100-x)"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    EXPECT_DOUBLE_EQ(joins[0].x, 100.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::Singular);
}

TEST(DetectPotentialStructureTest, HandlesInteriorJoinWithInfiniteWidthPiece)
{
    // Right piece [5, inf) has infinite width, exercising
    // isSingularApproaching's isfinite(pieceWidth)==false fallback (falls
    // back to the plain scale=max(|x0|,1) probe, matching A1's own formula,
    // since there's no finite piece boundary to clamp against). -1/x^2 is
    // smooth at x=5 (only singular at x=0, outside this piece), so the
    // join is a genuine Step (0 vs -1/25), not Singular -- confirming the
    // fallback doesn't spuriously flag a non-singular infinite-width piece.
    std::map<std::string, std::string> potential = {
        {"[0,5)", "0"}, {"[5, inf)", "-1/x^2"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    EXPECT_DOUBLE_EQ(joins[0].x, 5.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::Step);
}

// ---------------------------------------------------------------------------
// strategicKnotsFromJoins
// ---------------------------------------------------------------------------

TEST(StrategicKnotsFromJoinsTest, MultiplicityFormulaPerOrder)
{
    for (int order : {4, 6, 8})
    {
        std::vector<tise::DetectedJoin> joins = {
            {1.0, tise::JoinType::Step}, {2.0, tise::JoinType::StitchedKink}
        };
        auto knots = tise::strategicKnotsFromJoins(joins, order);

        auto findKnotAt = [&](double x) -> const tise::StrategicKnot * {
            for (const auto &k : knots)
                if (k.x == x) return &k;
            return nullptr;
        };

        const tise::StrategicKnot *step = findKnotAt(1.0);
        ASSERT_NE(step, nullptr);
        EXPECT_EQ(step->extraMultiplicity, order - 3);

        const tise::StrategicKnot *kink = findKnotAt(2.0);
        if (order - 4 > 0)
        {
            ASSERT_NE(kink, nullptr);
            EXPECT_EQ(kink->extraMultiplicity, order - 4);
        }
        else
        {
            EXPECT_EQ(kink, nullptr);
        }
    }
}

TEST(StrategicKnotsFromJoinsTest, StitchedKinkClampsToZeroAtOrderFour)
{
    std::vector<tise::DetectedJoin> joins = {{3.0, tise::JoinType::StitchedKink}};
    auto knots = tise::strategicKnotsFromJoins(joins, 4);
    EXPECT_TRUE(knots.empty());
}

TEST(StrategicKnotsFromJoinsTest, SingularAndContinuousProduceNoKnots)
{
    std::vector<tise::DetectedJoin> joins = {
        {1.0, tise::JoinType::Singular}, {2.0, tise::JoinType::Continuous}
    };
    auto knots = tise::strategicKnotsFromJoins(joins, 8);
    EXPECT_TRUE(knots.empty());
}

// ---------------------------------------------------------------------------
// buildStrategicRadialGrid
// ---------------------------------------------------------------------------

TEST(BuildStrategicRadialGridTest, ExistingPointGetsExtraMultiplicityOnly)
{
    // nNodes=11 on [0,10] -> spacing 1.0; x=5.0 already on the grid.
    std::vector<tise::StrategicKnot> knots = {{5.0, 3}};
    auto grid = tise::buildStrategicRadialGrid(11, 0.0, 10.0, knots);
    EXPECT_EQ(grid.size(), 11u + 3u);
    EXPECT_EQ(std::count(grid.begin(), grid.end(), 5.0), 1 + 3);
}

TEST(BuildStrategicRadialGridTest, NewPointInsertionSplicesAndAddsMultiplicity)
{
    // nNodes=11 on [0,10] -> spacing 1.0; x=5.5 is NOT on the grid.
    std::vector<tise::StrategicKnot> knots = {{5.5, 2}};
    auto grid = tise::buildStrategicRadialGrid(11, 0.0, 10.0, knots);
    EXPECT_EQ(grid.size(), 11u + 1u + 2u);
    EXPECT_EQ(std::count(grid.begin(), grid.end(), 5.5), 1 + 2);
    EXPECT_TRUE(std::is_sorted(grid.begin(), grid.end()));
}

TEST(BuildStrategicRadialGridTest, RemainsNonDecreasing)
{
    std::vector<tise::StrategicKnot> knots = {{2.0, 4}, {5.5, 3}, {8.0, 1}};
    auto grid = tise::buildStrategicRadialGrid(11, 0.0, 10.0, knots);
    EXPECT_TRUE(std::is_sorted(grid.begin(), grid.end()));
}

TEST(BuildStrategicRadialGridTest, ReturnedSizeAccountsForAllInsertions)
{
    // One on-grid knot (5.0, +3) and one off-grid knot (5.5, +2): total
    // extra entries = 3 (multiplicity only) + (1 splice + 2 multiplicity).
    std::vector<tise::StrategicKnot> knots = {{5.0, 3}, {5.5, 2}};
    auto grid = tise::buildStrategicRadialGrid(11, 0.0, 10.0, knots);
    EXPECT_EQ(grid.size(), 11u + 3u + (1u + 2u));
}

// ---------------------------------------------------------------------------
// Strategic node placement -- literal "Done when" criterion
// ---------------------------------------------------------------------------

TEST(StrategicNodePlacementAccuracyTest, ImprovesOverUniformGridForBoxBarrier)
{
    // Particle in a box [0,10] with a rectangular barrier on [5,6]
    // (docs/TDSE-original-design/2026-06-28-config-yaml-schema-design.md's
    // own worked example -- the task-breakdown's literal "Done when" case).
    // Ground-state eigenvalue only; L is unused by fillBandedMatrices' own
    // potential evaluation (kept only for eigenvalueError's hydrogen-
    // specific analytic comparison, not used here), passed as 0.
    std::map<std::string, std::string> potential = {
        {"[0,5)", "0"}, {"[5,6]", "10"}, {"(6,10]", "0"}
    };
    const int order = 8;
    const int L = 0;
    const double rMin = 0.0, rMax = 10.0;

    auto groundStateEnergy = [&](const std::vector<double> &grid, int nNodesForInit) {
        bspline::BSpline bs;
        int info = bs.init(nNodesForInit, order, grid);
        EXPECT_EQ(info, 0);
        int nEn = bs.getNBSplines() - 2;
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L, potential);
        auto result = tise::solveGeneralizedEigenproblem(std::move(H), std::move(S), nEn, order);
        return result.values[0];
    };

    // Fine uniform reference: converged ground truth.
    const int nFine = 201;
    double eFine = groundStateEnergy(tise::buildUniformRadialGrid(nFine, rMin, rMax), nFine);

    // Coarse uniform: 21 nodes, spacing 0.5 -- x=5 and x=6 already fall
    // exactly on the grid, deliberately, so the strategic comparison below
    // adds only knot degeneracy, no new distinct spatial points -- the
    // cleanest possible reading of the task's "same node count" wording.
    const int nCoarse = 21;
    double eCoarseUniform = groundStateEnergy(tise::buildUniformRadialGrid(nCoarse, rMin, rMax), nCoarse);

    // Coarse strategic: same 21-node uniform base, degenerate knots added
    // at the two Step joins (x=5, x=6).
    auto joins = tise::detectPotentialStructure(potential);
    auto knots = tise::strategicKnotsFromJoins(joins, order);
    auto strategicGrid = tise::buildStrategicRadialGrid(nCoarse, rMin, rMax, knots);
    // Per buildStrategicRadialGrid's contract: pass the returned grid's OWN
    // size, not nCoarse -- see Gap 1 above.
    double eCoarseStrategic = groundStateEnergy(strategicGrid, static_cast<int>(strategicGrid.size()));

    double errUniform = std::abs(eCoarseUniform - eFine);
    double errStrategic = std::abs(eCoarseStrategic - eFine);
    EXPECT_LT(errStrategic, errUniform);
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
// fillBandedMatrices -- generalized drop-set (A4b)
// ---------------------------------------------------------------------------

TEST(FillBandedMatricesDropSetTest, ExplicitClassicDropSetMatchesDefault)
{
    const int order = 8, nNodes = 15;
    std::vector<double> grid(nNodes);
    for (int i = 0; i < nNodes; ++i) grid[i] = i;
    bspline::BSpline bs;
    ASSERT_EQ(bs.init(nNodes, order, grid), 0);
    const int nBSplines = bs.getNBSplines();
    const int nEn = nBSplines - 2;

    std::map<std::string, std::string> potential = {{"[0," + std::to_string(nNodes - 1) + "]", "x"}};
    auto [H1, S1] = tise::fillBandedMatrices(bs, nEn, order, 0, potential);
    auto [H2, S2] = tise::fillBandedMatrices(bs, nEn, order, 0, potential, std::vector<int>{1, nBSplines});

    EXPECT_EQ(H1, H2);
    EXPECT_EQ(S1, S2);
}

TEST(FillBandedMatricesDropSetTest, InteriorDropSetMatchesDirectIntegralAtMappedPositions)
{
    // order=4, nNodes=6 -> nBSplines=8. dropSet={1,4,8}: both classic ends
    // plus one interior index. Kept={2,3,5,6,7} (nEn=5) get columns
    // 1,2,3,4,5 respectively (hand-derived: docs/planning/engineer-a-plan-A4b.md).
    const int order = 4, nNodes = 6;
    std::vector<double> grid = {0, 1, 2, 3, 4, 5};
    bspline::BSpline bs;
    ASSERT_EQ(bs.init(nNodes, order, grid), 0);
    ASSERT_EQ(bs.getNBSplines(), 8);

    std::map<std::string, std::string> potential = {{"[0,5]", "1.0"}};
    std::vector<int> dropSet = {1, 4, 8};
    const int nEn = 5;
    auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, 0, potential, dropSet);

    bspline::D2DFun fUni = [](double, const double *) { return 1.0; };
    bspline::D2DFun fPot = [&](double x, const double *) { return tise::evaluateFunction(potential, x); };
    double parvec[1] = {0.0};
    auto bandIndex = [&](int row, int col) { return (row - 1) + (col - 1) * order; };
    auto directH = [&](int iBs1, int iBs2) {
        return bs.integral(fUni, iBs1, iBs2, 1, 1) / 2.0 + bs.integral(fPot, iBs1, iBs2, 0, 0, parvec);
    };

    // Physical (5,5) -> col(3,3), row=4+3-3=4, idx=bandIndex(4,3)=11.
    EXPECT_NEAR(H[bandIndex(4, 3)], directH(5, 5), 1e-12);
    // Physical (3,5) -> col(2,3), row=4+2-3=3, idx=bandIndex(3,3)=10.
    EXPECT_NEAR(H[bandIndex(3, 3)], directH(3, 5), 1e-12);
    // Physical (2,5) -> col(1,3), row=4+1-3=2, idx=bandIndex(2,3)=9.
    EXPECT_NEAR(H[bandIndex(2, 3)], directH(2, 5), 1e-12);
}

TEST(FillBandedMatricesDropSetTest, EmptyDropSetKeepsAllBSplines)
{
    const int order = 4, nNodes = 6;
    std::vector<double> grid = {0, 1, 2, 3, 4, 5};
    bspline::BSpline bs;
    ASSERT_EQ(bs.init(nNodes, order, grid), 0);
    const int nBSplines = bs.getNBSplines();
    ASSERT_EQ(nBSplines, 8);

    std::map<std::string, std::string> potential = {{"[0,5]", "1.0"}};
    auto [H, S] = tise::fillBandedMatrices(bs, nBSplines, order, 0, potential, std::vector<int>{});

    // Physical B-spline 1 (normally dropped) now has column 1; its diagonal
    // overlap S(1,1) = ||B_1||^2 > 0 must appear at bandIndex(order,1).
    auto bandIndex = [&](int row, int col) { return (row - 1) + (col - 1) * order; };
    EXPECT_GT(S[bandIndex(order, 1)], 0.0);
}

TEST(FillBandedMatricesDropSetTest, MismatchedNEnThrows)
{
    const int order = 4, nNodes = 6;
    std::vector<double> grid = {0, 1, 2, 3, 4, 5};
    bspline::BSpline bs;
    ASSERT_EQ(bs.init(nNodes, order, grid), 0);
    std::map<std::string, std::string> potential = {{"[0,5]", "1.0"}};
    // nBSplines=8, dropSet={1,8} -> true nEn=6, but pass 5.
    EXPECT_THROW(
        tise::fillBandedMatrices(bs, 5, order, 0, potential, std::vector<int>{1, 8}),
        std::runtime_error);
}

TEST(FillBandedMatricesDropSetTest, OutOfRangeDropSetIndexThrows)
{
    const int order = 4, nNodes = 6;
    std::vector<double> grid = {0, 1, 2, 3, 4, 5};
    bspline::BSpline bs;
    ASSERT_EQ(bs.init(nNodes, order, grid), 0);
    std::map<std::string, std::string> potential = {{"[0,5]", "1.0"}};
    EXPECT_THROW(
        tise::fillBandedMatrices(bs, 7, order, 0, potential, std::vector<int>{0, 8}),
        std::runtime_error);
}

// ---------------------------------------------------------------------------
// eigenstateCoefficients -- generalized drop-set (A4b)
// ---------------------------------------------------------------------------

TEST(EigenstateCoefficientsDropSetTest, InteriorDropSetZeroesExactlyDroppedIndices)
{
    // Same nBSplines=8/dropSet={1,4,8} as the fillBandedMatrices test above.
    // Synthetic single-column evec: components 10,20,30,40,50 for the 5
    // kept columns.
    const int nBSplines = 8, nEn = 5;
    std::vector<double> evec = {10, 20, 30, 40, 50};
    std::vector<int> dropSet = {1, 4, 8};

    auto coeffs = tise::eigenstateCoefficients(evec, 1, nEn, nBSplines, dropSet);
    ASSERT_EQ(coeffs.size(), 8u);

    EXPECT_DOUBLE_EQ(coeffs[0], 0.0); // physical 1, dropped
    EXPECT_DOUBLE_EQ(coeffs[3], 0.0); // physical 4, dropped
    EXPECT_DOUBLE_EQ(coeffs[7], 0.0); // physical 8, dropped
    EXPECT_DOUBLE_EQ(coeffs[1], 10.0); // physical 2 -> column 1
    EXPECT_DOUBLE_EQ(coeffs[2], 20.0); // physical 3 -> column 2
    EXPECT_DOUBLE_EQ(coeffs[4], 30.0); // physical 5 -> column 3
    EXPECT_DOUBLE_EQ(coeffs[5], 40.0); // physical 6 -> column 4
    EXPECT_DOUBLE_EQ(coeffs[6], 50.0); // physical 7 -> column 5
}

TEST(EigenstateCoefficientsDropSetTest, ClassicDropSetMatchesDefault)
{
    const int nBSplines = 8, nEn = 6;
    std::vector<double> evec = {1, 2, 3, 4, 5, 6};
    auto c1 = tise::eigenstateCoefficients(evec, 1, nEn, nBSplines);
    auto c2 = tise::eigenstateCoefficients(evec, 1, nEn, nBSplines, std::vector<int>{1, nBSplines});
    EXPECT_EQ(c1, c2);
}

TEST(EigenstateCoefficientsDropSetTest, MismatchedNEnThrows)
{
    std::vector<double> evec(5, 0.0);
    EXPECT_THROW(
        tise::eigenstateCoefficients(evec, 1, 5, 8, std::vector<int>{1, 8}), // true nEn=6
        std::runtime_error);
}

// ---------------------------------------------------------------------------
// bSplinesTouchingX
// ---------------------------------------------------------------------------

TEST(BSplinesTouchingXTest, MatchesHandComputedSupportAtInteriorPoint)
{
    std::vector<double> grid(11);
    for (int i = 0; i < 11; ++i) grid[i] = i;
    auto touching = tise::bSplinesTouchingX(11, 4, grid, 4.5);
    std::vector<int> expected = {5, 6, 7, 8};
    EXPECT_EQ(touching, expected);
}

TEST(BSplinesTouchingXTest, MatchesHandComputedSupportAtGridBoundary)
{
    std::vector<double> grid(11);
    for (int i = 0; i < 11; ++i) grid[i] = i;
    auto touching = tise::bSplinesTouchingX(11, 4, grid, 5.0);
    std::vector<int> expected = {5, 6, 7, 8, 9};
    EXPECT_EQ(touching, expected);
}

TEST(BSplinesTouchingXTest, EmptyOutsideDomain)
{
    std::vector<double> grid(11);
    for (int i = 0; i < 11; ++i) grid[i] = i;
    auto touching = tise::bSplinesTouchingX(11, 4, grid, 15.0);
    EXPECT_TRUE(touching.empty());
}

// ---------------------------------------------------------------------------
// Round-trip and end-to-end (A4b's actual purpose: detection -> removal)
// ---------------------------------------------------------------------------

TEST(DropSetRoundTripTest, CoefficientsZeroExactlyAtDroppedIndicesAfterSolve)
{
    const int order = 8, nNodes = 21;
    std::vector<double> grid(nNodes);
    for (int i = 0; i < nNodes; ++i) grid[i] = i * 10.0 / (nNodes - 1);
    bspline::BSpline bs;
    ASSERT_EQ(bs.init(nNodes, order, grid), 0);
    const int nBSplines = bs.getNBSplines();

    std::map<std::string, std::string> potential = {{"[0,10]", "0"}};
    const int interior = nBSplines / 2;
    std::vector<int> dropSet = {1, interior, nBSplines};
    const int nEn = nBSplines - 3;

    auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, 0, potential, dropSet);
    auto result = tise::solveGeneralizedEigenproblem(std::move(H), std::move(S), nEn, order);
    ASSERT_EQ(result.dim, nEn);

    auto coeffs = tise::eigenstateCoefficients(result.vectors, 1, nEn, nBSplines, dropSet);
    EXPECT_DOUBLE_EQ(coeffs[0], 0.0);
    EXPECT_DOUBLE_EQ(coeffs[interior - 1], 0.0);
    EXPECT_DOUBLE_EQ(coeffs[nBSplines - 1], 0.0);
}

TEST(SingularPotentialBSplineRemovalTest, RemovingBSplinesNearSingularityProducesFiniteSolution)
{
    // Coulomb-like singular potential at x=0 (mirrors A4a's own
    // DetectPotentialStructureTest.SingularAtCoulombOrigin). Runs the full
    // detection-to-removal pipeline: detectPotentialStructure finds the
    // Singular join, bSplinesTouchingX turns it into a drop-set, the
    // generalized fillBandedMatrices/eigenstateCoefficients use it -- this
    // is what makes A4a's Singular detection "more than detection-only"
    // per the master plan's own framing of A4b's purpose.
    const int order = 8, nNodes = 41;
    const double rMin = 0.0, rMax = 40.0;
    std::map<std::string, std::string> potential = {{"(0, 40]", "-1/x + 1/x^2"}};

    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    ASSERT_EQ(joins[0].type, tise::JoinType::Singular);
    ASSERT_DOUBLE_EQ(joins[0].x, 0.0);

    auto grid = tise::buildUniformRadialGrid(nNodes, rMin, rMax);
    // x=0 is the domain's own left boundary, not an interior point -- pulls
    // in the entire first cluster of `order` B-splines (hand-traced above),
    // not just one.
    auto dropSet = tise::bSplinesTouchingX(nNodes, order, grid, joins[0].x);
    ASSERT_EQ(dropSet.size(), static_cast<size_t>(order));

    bspline::BSpline bs;
    ASSERT_EQ(bs.init(nNodes, order, grid), 0);
    const int nBSplines = bs.getNBSplines();
    const int nEn = nBSplines - static_cast<int>(dropSet.size());

    auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, 0, potential, dropSet);
    auto result = tise::solveGeneralizedEigenproblem(std::move(H), std::move(S), nEn, order);

    ASSERT_EQ(result.values.size(), static_cast<size_t>(nEn));
    for (double e : result.values)
        EXPECT_TRUE(std::isfinite(e));

    auto coeffs = tise::eigenstateCoefficients(result.vectors, 1, nEn, nBSplines, dropSet);
    for (int bsIdx : dropSet)
        EXPECT_DOUBLE_EQ(coeffs[bsIdx - 1], 0.0);
}

// ---------------------------------------------------------------------------
// computeEAcc
// ---------------------------------------------------------------------------

TEST(ComputeEAccTest, MatchesClosedForm)
{
    EXPECT_NEAR(tise::computeEAcc(0.5, 1.0), 2.0 * M_PI * M_PI, 1e-12);
}

TEST(ComputeEAccTest, ScalesInverselyWithMassAndSpacingSquared)
{
    double base = tise::computeEAcc(0.5, 1.0);

    // Doubling nodeSpacing (mass fixed): E_acc ~ 1/dx^2 -> quarters.
    EXPECT_NEAR(tise::computeEAcc(1.0, 1.0), base / 4.0, 1e-12);

    // Doubling mass (nodeSpacing fixed): E_acc ~ 1/m -> halves.
    EXPECT_NEAR(tise::computeEAcc(0.5, 2.0), base / 2.0, 1e-12);
}

// ---------------------------------------------------------------------------
// warnIfContinuumExceedsEAcc
// ---------------------------------------------------------------------------

TEST(WarnIfContinuumExceedsEAccTest, FiresWhenEMaxExceeds)
{
    std::ostringstream warn;
    EXPECT_TRUE(tise::warnIfContinuumExceedsEAcc(10.0, 5.0, warn));
    EXPECT_FALSE(warn.str().empty());
}

TEST(WarnIfContinuumExceedsEAccTest, SilentWhenEMaxBelow)
{
    std::ostringstream warn;
    EXPECT_FALSE(tise::warnIfContinuumExceedsEAcc(3.0, 5.0, warn));
    EXPECT_TRUE(warn.str().empty());
}

TEST(WarnIfContinuumExceedsEAccTest, BoundaryEqualDoesNotFire)
{
    std::ostringstream warn;
    EXPECT_FALSE(tise::warnIfContinuumExceedsEAcc(5.0, 5.0, warn));
    EXPECT_TRUE(warn.str().empty());
}

TEST(WarnIfContinuumExceedsEAccTest, WarningTextMentionsUnreliableResults)
{
    std::ostringstream warn;
    tise::warnIfContinuumExceedsEAcc(10.0, 5.0, warn);
    EXPECT_NE(warn.str().find("unreliable"), std::string::npos);
}
