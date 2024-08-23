// main.cpp
// This is an open source implementation of the Schumaker splines in C++. Copyright Stuart Baumann under the MIT License (https://en.wikipedia.org/wiki/MIT_License).
// The algorithmic logic of the spline is from "Schumaker, L.L. 1983. On shape-preserving quadratic spline interpolation. SIAM Journal of Numerical Analysis 20: 854-64." and is also described in the book "Numerical Methods in Economics" by Kenneth L. Judd.
// For an explanation of why shape preserving splines are useful for economic dynamics problems see https://cran.r-project.org/web/packages/schumaker/vignettes/schumaker.html.
// The original code was written by S. Baumann and can be found at https://github.com/s-baumann/Schumaker.cpp. Pull requests welcome.

// Note that this file just demonstrates the usage of the Schumaker spline. The actual implementation is in Schumaker.cpp and Schumaker.hpp.

#include <cmath>
#include <iostream>
#include <vector>

#include "Schumaker.hpp"
#include "Schumaker.cpp"
template <class T> std::string vec2string(const std::vector<T> v) {
  std::string s = "std::vector{T} {";
  for (int i = 0; i < v.size(); i++) {
    if (i > 0) {
      s = s + ", ";
    }
    s = s + std::to_string(v[i]);
  }
  s = s + " }";
  return s;
}

int main() {
  // Making test data
  const auto x_vec = std::vector<double>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<double> y(x_vec.size());
  for (auto i = 0; i < x_vec.size(); i++) {
    y[i] = (std::pow(x_vec[i], 2));
  }

  // Constructing a spline
  auto spl = Spline(x_vec, y);
  // Evaluating the spline and its derivatives.
  const auto testx =
      std::vector<double>{0, 0.5, 1, 2, 3, 4, 5, 6, 6.6, 7, 8, 9, 10, 12};
  for (auto x : testx) {
    std::cout << " at input of " << x << " we have " << spl.spline(x)
              << " for the value, " << spl.derivative(x)
              << " for the first derivative, " << spl.second_derivative(x)
              << " for the second derivative, " << std::endl;
  }

  // Inputting a vector to the spline.
  const auto result0 = spl.spline(testx);
  const auto result1 = spl.derivative(testx);
  const auto result2 = spl.second_derivative(testx);
  std::cout << " at input of " << vec2string(testx) << " we have " << vec2string(result0)
              << " for the value, " << vec2string(result1)
              << " for the first derivative, " << vec2string(result2)
              << " for the second derivative, " << std::endl;

  return 0;
}