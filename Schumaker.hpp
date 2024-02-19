
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

  double
  spline(double inp);

  double derivative(
      double inp);

  double second_derivative(double inp);

  std::vector<double> Intervals_;
  std::vector<std::vector<double>> Coefficients_;
};
