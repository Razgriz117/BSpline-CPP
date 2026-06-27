#include <vector>
#include <string>
#include "BSpline.hpp"

void printBandMatrixAsDense(
     const std::vector<bspline::Real> &band,
     int n,
     int kd,
     int ld,
     const std::string &name,
     int width = 14,
     int precision = 6);

std::vector<bspline::Real> buildAdjointMatrix(
     const std::vector<bspline::Real> &evec,
     int n,
     int ldz
);