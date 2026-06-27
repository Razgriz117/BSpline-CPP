#include <gtest/gtest.h>
#include "utils.hpp"

#include <vector>
#include <sstream>
#include <iostream>

using Real = bspline::Real;

// ---------------------------------------------------------------------------
// buildAdjointMatrix tests
// ---------------------------------------------------------------------------

static std::vector<Real> identityColMajor(int n)
{
    std::vector<Real> m(n * n, 0.0);
    for (int i = 0; i < n; ++i)
        m[i + i * n] = 1.0;
    return m;
}

TEST(BuildAdjointMatrixTest, IdentityTransposeIsIdentity)
{
    const int n = 4;
    auto I = identityColMajor(n);
    auto C = buildAdjointMatrix(I, n, n);

    ASSERT_EQ(static_cast<int>(C.size()), n * n);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
        {
            Real expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(C[i + j * n], expected, 1e-15)
                << "C[" << i << "," << j << "] wrong";
        }
}

TEST(BuildAdjointMatrixTest, KnownValues2x2)
{
    // Z = [[1,2],[3,4]] stored column-major: {1,3, 2,4}
    // Z^T = [[1,3],[2,4]] stored column-major: {1,2, 3,4}
    int n = 2, ldz = 2;
    std::vector<Real> Z = {1.0, 3.0,   // column 0
                            2.0, 4.0};  // column 1
    auto C = buildAdjointMatrix(Z, n, ldz);

    ASSERT_EQ(static_cast<int>(C.size()), 4);
    EXPECT_NEAR(C[0 + 0*2], 1.0, 1e-15); // C(0,0) = Z(0,0) = 1
    EXPECT_NEAR(C[1 + 0*2], 2.0, 1e-15); // C(1,0) = Z(0,1) = 2
    EXPECT_NEAR(C[0 + 1*2], 3.0, 1e-15); // C(0,1) = Z(1,0) = 3
    EXPECT_NEAR(C[1 + 1*2], 4.0, 1e-15); // C(1,1) = Z(1,1) = 4
}

TEST(BuildAdjointMatrixTest, KnownValues3x3)
{
    // Z[i + j*ldz] = i + j*ldz (0-based)
    int n = 3, ldz = 3;
    std::vector<Real> Z(9);
    for (int k = 0; k < 9; ++k)
        Z[k] = static_cast<Real>(k);

    auto C = buildAdjointMatrix(Z, n, ldz);

    // C[i + j*n] = Z[j + i*ldz]
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            EXPECT_NEAR(C[i + j * n],
                        static_cast<Real>(j + i * ldz),
                        1e-15)
                << "C[" << i << "," << j << "] wrong";
}

TEST(BuildAdjointMatrixTest, PaddedLdz_IgnoresExtraRows)
{
    // Z is 4×2 (ldz=4) but only top 2×2 is used (n=2)
    int n = 2, ldz = 4;
    std::vector<Real> Z = {1.0, 2.0, 99.0, 99.0,   // column 0
                            3.0, 4.0, 99.0, 99.0};  // column 1
    auto C = buildAdjointMatrix(Z, n, ldz);

    ASSERT_EQ(static_cast<int>(C.size()), 4);
    EXPECT_NEAR(C[0 + 0*2], 1.0, 1e-15); // C(0,0) = Z[0 + 0*4] = 1
    EXPECT_NEAR(C[1 + 0*2], 3.0, 1e-15); // C(1,0) = Z[0 + 1*4] = 3
    EXPECT_NEAR(C[0 + 1*2], 2.0, 1e-15); // C(0,1) = Z[1 + 0*4] = 2
    EXPECT_NEAR(C[1 + 1*2], 4.0, 1e-15); // C(1,1) = Z[1 + 1*4] = 4
    for (Real v : C)
        EXPECT_NE(v, 99.0) << "Padding value leaked into result";
}

TEST(BuildAdjointMatrixTest, LargerMatrix)
{
    int n = 5, ldz = 5;
    std::vector<Real> Z(25);
    for (int k = 0; k < 25; ++k)
        Z[k] = static_cast<Real>(k + 1);

    auto C = buildAdjointMatrix(Z, n, ldz);

    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            EXPECT_NEAR(C[i + j * n], Z[j + i * ldz], 1e-15)
                << "Mismatch at (" << i << "," << j << ")";
}

// ---------------------------------------------------------------------------
// printBandMatrixAsDense smoke tests (verify no crash, suppress output)
// ---------------------------------------------------------------------------

TEST(PrintBandMatrixTest, SmokePrint2x2Band)
{
    // 2×2 symmetric band, kd=1, ld=2
    int n = 2, kd = 1, ld = 2;
    std::vector<Real> band(ld * n, 0.0);
    // A(1,1) at AB(kd+1,1): row=kd+1=2, col=1 => idx=(2-1)+(1-1)*ld=1
    band[1] = 1.0;
    // A(1,2) at AB(kd+1-1,2)=AB(1,2): row=1, col=2 => idx=(1-1)+(2-1)*ld=2
    band[2] = 2.0;
    // A(2,2) at AB(kd+1,2)=AB(2,2): row=2, col=2 => idx=(2-1)+(2-1)*ld=3
    band[3] = 3.0;

    std::streambuf *old = std::cout.rdbuf();
    std::ostringstream sink;
    std::cout.rdbuf(sink.rdbuf());

    EXPECT_NO_THROW(printBandMatrixAsDense(band, n, kd, ld, "Test2x2"));

    std::cout.rdbuf(old);
}

TEST(PrintBandMatrixTest, SmokePrint4x4Band)
{
    int n = 4, kd = 2, ld = 3;
    std::vector<Real> band(ld * n, 0.0);
    // Fill the main diagonal: A(j,j) stored at row=kd+1, col=j
    // idx = (kd+1-1) + (j-1)*ld = kd + (j-1)*ld
    for (int j = 1; j <= n; ++j)
        band[kd + (j - 1) * ld] = static_cast<Real>(j);

    std::streambuf *old = std::cout.rdbuf();
    std::ostringstream sink;
    std::cout.rdbuf(sink.rdbuf());

    EXPECT_NO_THROW(printBandMatrixAsDense(band, n, kd, ld, "Diag4x4"));

    std::cout.rdbuf(old);
}
