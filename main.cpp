#include <cmath>
#include <iostream>
#include <vector>

#include "Schumaker.hpp"

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
  const auto &x = std::vector<double>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<double> y(x.size());
  for (auto i = 0; i < x.size(); i++) {
    y[i] = (std::pow(x[i], 2));
  }

  // Constructing a spline
  auto spl = Spline(x, y);

  // Evaluating the spline and its derivatives.
  const auto &testx =
      std::vector<double>{0, 0.5, 1, 2, 3, 4, 5, 6, 6.6, 7, 8, 9, 10, 12};
  for (const auto &x : testx) {
    std::cout << " at input of " << x << " we have " << spl.spline(x)
              << " for the value, " << spl.derivative(x)
              << " for the first derivative, " << spl.second_derivative(x)
              << " for the second derivative, " << std::endl;
  }

  return 0;
}