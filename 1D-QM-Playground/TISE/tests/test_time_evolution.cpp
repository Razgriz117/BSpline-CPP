#include <gtest/gtest.h>
#include "time_evolution.hpp"
#include "tise.hpp"
#include "BSpline.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <sstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// gaussianWavepacket
// ---------------------------------------------------------------------------

TEST(GaussianWavepacketTest, PeakIsAtR0)
{
    const double r0 = 5.0, mOmega = 0.5, hbar = 1.0;
    double peak = tevol::gaussianWavepacket(r0, r0, mOmega, hbar);
    // Values either side should be smaller
    EXPECT_GT(peak, tevol::gaussianWavepacket(r0 - 1.0, r0, mOmega, hbar));
    EXPECT_GT(peak, tevol::gaussianWavepacket(r0 + 1.0, r0, mOmega, hbar));
}

TEST(GaussianWavepacketTest, SymmetricAroundR0)
{
    const double r0 = 5.0, mOmega = 0.5, hbar = 1.0, dx = 0.7;
    double left  = tevol::gaussianWavepacket(r0 - dx, r0, mOmega, hbar);
    double right = tevol::gaussianWavepacket(r0 + dx, r0, mOmega, hbar);
    EXPECT_NEAR(left, right, 1e-15);
}

TEST(GaussianWavepacketTest, PeakValueIsPositive)
{
    EXPECT_GT(tevol::gaussianWavepacket(3.0, 3.0, 0.5, 1.0), 0.0);
}

TEST(GaussianWavepacketTest, DecaysAwayFromPeak)
{
    const double r0 = 5.0, mOmega = 1.0, hbar = 1.0;
    double v0 = tevol::gaussianWavepacket(r0,       r0, mOmega, hbar);
    double v1 = tevol::gaussianWavepacket(r0 + 2.0, r0, mOmega, hbar);
    double v2 = tevol::gaussianWavepacket(r0 + 4.0, r0, mOmega, hbar);
    EXPECT_GT(v0, v1);
    EXPECT_GT(v1, v2);
}

// ---------------------------------------------------------------------------
// timeEvolveState
// ---------------------------------------------------------------------------

class TimeEvolveStateTest : public ::testing::Test
{
protected:
    int n = 4;
    std::vector<double> eval = {-0.5, -0.125, -1.0/18.0, -1.0/32.0};
    Eigen::VectorXd Phi_G;
    double hbar = 1.0;

    void SetUp() override
    {
        Phi_G = Eigen::VectorXd::Ones(n) * 2.0;
    }
};

TEST_F(TimeEvolveStateTest, AtTZeroEqualsPhiG)
{
    auto V = tevol::timeEvolveState(Phi_G, eval, 0.0, hbar);
    ASSERT_EQ(V.size(), n);
    for (int i = 0; i < n; ++i)
    {
        EXPECT_NEAR(V(i).real(), Phi_G(i), 1e-14);
        EXPECT_NEAR(V(i).imag(), 0.0,       1e-14);
    }
}

TEST_F(TimeEvolveStateTest, MagnitudesPreservedAtLaterTime)
{
    auto V = tevol::timeEvolveState(Phi_G, eval, 1.5, hbar);
    for (int i = 0; i < n; ++i)
        EXPECT_NEAR(std::abs(V(i)), std::abs(Phi_G(i)), 1e-13);
}

TEST_F(TimeEvolveStateTest, PhaseMatchesExplicitFormula)
{
    double t = 0.3;
    auto V = tevol::timeEvolveState(Phi_G, eval, t, hbar);
    for (int i = 0; i < n; ++i)
    {
        std::complex<double> expected =
            std::exp(std::complex<double>(0.0, -eval[i] * t / hbar)) * Phi_G(i);
        EXPECT_NEAR(V(i).real(), expected.real(), 1e-13);
        EXPECT_NEAR(V(i).imag(), expected.imag(), 1e-13);
    }
}

// ---------------------------------------------------------------------------
// projectToEigenBasis
// ---------------------------------------------------------------------------

TEST(ProjectToEigenBasisTest, IdentityMatrixReturnsBG)
{
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);
    Eigen::VectorXd B_G(3);
    B_G << 1.0, 2.0, 3.0;
    auto Phi = tevol::projectToEigenBasis(I, B_G);
    for (int i = 0; i < 3; ++i)
        EXPECT_NEAR(Phi(i), B_G(i), 1e-15);
}

TEST(ProjectToEigenBasisTest, Known2x2)
{
    // Cinv = [[2,0],[0,3]], B_G = [1,1] => Phi = [2,3]
    Eigen::MatrixXd Cinv(2, 2);
    Cinv << 2.0, 0.0, 0.0, 3.0;
    Eigen::VectorXd B_G(2);
    B_G << 1.0, 1.0;
    auto Phi = tevol::projectToEigenBasis(Cinv, B_G);
    EXPECT_NEAR(Phi(0), 2.0, 1e-15);
    EXPECT_NEAR(Phi(1), 3.0, 1e-15);
}

// ---------------------------------------------------------------------------
// transformToSpaceBasis
// ---------------------------------------------------------------------------

TEST(TransformToSpaceBasisTest, IdentityMatrixPreservesVt)
{
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);
    Eigen::VectorXcd V_t(3);
    V_t << std::complex<double>(1, 2), std::complex<double>(3, 4), std::complex<double>(5, 6);
    auto CV = tevol::transformToSpaceBasis(I, V_t);
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(CV(i).real(), V_t(i).real(), 1e-14);
        EXPECT_NEAR(CV(i).imag(), V_t(i).imag(), 1e-14);
    }
}

TEST(TransformToSpaceBasisTest, Known2x2)
{
    // C = [[1,2],[3,4]], V_t = [1+0i, 0+0i] => CV = [1, 3]
    Eigen::MatrixXd C(2, 2);
    C << 1.0, 2.0, 3.0, 4.0;
    Eigen::VectorXcd V_t(2);
    V_t << std::complex<double>(1, 0), std::complex<double>(0, 0);
    auto CV = tevol::transformToSpaceBasis(C, V_t);
    EXPECT_NEAR(CV(0).real(), 1.0, 1e-14);
    EXPECT_NEAR(CV(1).real(), 3.0, 1e-14);
}

// ---------------------------------------------------------------------------
// computeGaussianOverlaps
// ---------------------------------------------------------------------------

class ComputeGaussianOverlapsTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> grid(nNodes);
        for (int i = 0; i < nNodes; ++i)
            grid[i] = rMin + (rMax - rMin) * i / (nNodes - 1);
        ASSERT_EQ(bs.init(nNodes, order, grid), 0);
        int nBs = bs.getNBSplines();
        nEn = nBs - 2;
    }

    bspline::BSpline bs;
    int nNodes = 21, order = 6;
    double rMin = 0.0, rMax = 30.0;
    double r0 = 5.0, mOmega = 0.6847, hbar = 1.0;
    int nEn = 0;
};

TEST_F(ComputeGaussianOverlapsTest, LengthIsNEn)
{
    auto B_G = tevol::computeGaussianOverlaps(bs, nEn, order, r0, mOmega, hbar);
    EXPECT_EQ(B_G.size(), nEn);
}

TEST_F(ComputeGaussianOverlapsTest, AllFinite)
{
    auto B_G = tevol::computeGaussianOverlaps(bs, nEn, order, r0, mOmega, hbar);
    for (int i = 0; i < nEn; ++i)
        EXPECT_TRUE(std::isfinite(B_G(i)));
}

TEST_F(ComputeGaussianOverlapsTest, NonNegativeNearPeak)
{
    // B-splines near the Gaussian center should give positive overlap
    auto B_G = tevol::computeGaussianOverlaps(bs, nEn, order, r0, mOmega, hbar);
    // At least one element should be positive
    bool anyPositive = false;
    for (int i = 0; i < nEn; ++i)
        if (B_G(i) > 0.0) anyPositive = true;
    EXPECT_TRUE(anyPositive);
}

// ---------------------------------------------------------------------------
// writeTimestep
// ---------------------------------------------------------------------------

class WriteTimestepTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> grid(nNodes);
        for (int i = 0; i < nNodes; ++i)
            grid[i] = i * 1.0 / (nNodes - 1);
        ASSERT_EQ(bs.init(nNodes, order, grid), 0);
        nBSplines = bs.getNBSplines();
        CV_padded = Eigen::VectorXcd::Zero(nBSplines);
        CV_padded(nBSplines / 2) = std::complex<double>(1.0, 0.5);
    }

    bspline::BSpline bs;
    int nNodes = 11, order = 4, nBSplines = 0;
    Eigen::VectorXcd CV_padded;
};

TEST_F(WriteTimestepTest, LineCountEqualsNpts)
{
    std::ostringstream ss;
    tevol::writeTimestep(ss, bs, CV_padded, nBSplines, 20, 0.0, 1.0);
    std::string output = ss.str();
    int lines = static_cast<int>(std::count(output.begin(), output.end(), '\n'));
    EXPECT_EQ(lines, 20);
}

TEST_F(WriteTimestepTest, ThreeColumnsPerLine)
{
    std::ostringstream ss;
    tevol::writeTimestep(ss, bs, CV_padded, nBSplines, 5, 0.0, 1.0);
    std::istringstream in(ss.str());
    double a, b, c;
    int count = 0;
    while (in >> a >> b >> c) ++count;
    EXPECT_EQ(count, 5);
}

TEST_F(WriteTimestepTest, ValuesAreFinite)
{
    std::ostringstream ss;
    tevol::writeTimestep(ss, bs, CV_padded, nBSplines, 10, 0.0, 1.0);
    std::istringstream in(ss.str());
    double a, b, c;
    while (in >> a >> b >> c)
    {
        EXPECT_TRUE(std::isfinite(a));
        EXPECT_TRUE(std::isfinite(b));
        EXPECT_TRUE(std::isfinite(c));
    }
}
