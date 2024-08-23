// Schumaker.hpp
// This is an open source implementation of the Schumaker splines in C++. Copyright Stuart Baumann under the MIT License (https://en.wikipedia.org/wiki/MIT_License).
// The algorithmic logic of the spline is from "Schumaker, L.L. 1983. On shape-preserving quadratic spline interpolation. SIAM Journal of Numerical Analysis 20: 854-64." and is also described in the book "Numerical Methods in Economics" by Kenneth L. Judd.
// For an explanation of why shape preserving splines are useful for economic dynamics problems see https://cran.r-project.org/web/packages/schumaker/vignettes/schumaker.html.
// The original code was written by S. Baumann and can be found at https://github.com/s-baumann/Schumaker.cpp. Pull requests welcome.

#pragma once

#include <cfloat>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>


class Spline {
public:
  Spline(const std::vector<double> x,
         const std::vector<double> y);

  double spline(double inp);

  double derivative(double inp);

  double second_derivative(double inp);

  std::vector<double> spline(std::vector<double> inp);

  std::vector<double> derivative(std::vector<double> inp);

  std::vector<double> second_derivative(std::vector<double> inp);

  std::vector<double> Intervals_;
  std::vector<std::vector<double>> Coefficients_;
};
