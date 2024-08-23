// Schumaker.cpp
// This is an open source implementation of the Schumaker splines in C++. Copyright Stuart Baumann under the MIT License (https://en.wikipedia.org/wiki/MIT_License).
// The algorithmic logic of the spline is from "Schumaker, L.L. 1983. On shape-preserving quadratic spline interpolation. SIAM Journal of Numerical Analysis 20: 854-64." and is also described in the book "Numerical Methods in Economics" by Kenneth L. Judd.
// For an explanation of why shape preserving splines are useful for economic dynamics problems see https://cran.r-project.org/web/packages/schumaker/vignettes/schumaker.html.
// The original code was written by S. Baumann and can be found at https://github.com/s-baumann/Schumaker.cpp. Pull requests welcome.
#include "Schumaker.hpp"

#include <cfloat>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

double sum_vector(const std::vector<double> x) {
  return accumulate(x.begin(), x.end(), 0.0);
}

template <class T>
std::vector<T> append_vectors(const std::vector<T> v1,
                              const std::vector<T> v2) {
  std::vector<T> newvec = v1; // possible optimisation here if willing to make
                              // v1 nonconst and the code nonfunctional.
  for (const auto &i : v2) {
    newvec.push_back(i);
  }
  return newvec;
}

int search_sorted_last_less_than_x(const std::vector<double> vec,
                                   const double x) {
  for (int i = vec.size() - 1; i > 0; i--) {
    if (x >= vec[i]) {
      return i;
    }
  }
  return 0;
}

/*
    Splits an interval into 2 subintervals and
        creates the quadratic coefficients## #Inputs * `gradients` -
    A 2 entry vector with gradients at either end of the interval * `y` -
    A 2 entry vector with y values at either end of the interval * `x` -
    A 2 entry vector with x values at either end of the interval## #Returns *
        A 2 x 4 matrix
            .The first column is the x values of start of the two subintervals
            .The last 3 columns
                are quadratic coefficients in two subintervals.
*/

std::tuple<std::vector<double>, std::vector<std::vector<double>>>
schumakerIndInterval(const std::vector<double> gradients,
                     const std::vector<double> y, const std::vector<double> x) {
  double tsi;

  // The SchumakerIndInterval function takes in each interval individually
  // and returns the location of the knot as well as the quadratic
  // coefficients in \
//   each subinterval.
  // #Judd(1998), page 232, Lemma 6.11.1 provides this if condition:
  if (sum_vector(gradients) * (x[1] - x[0]) == 2 * (y[1] - y[0])) {
    tsi = x[1];
  } else {
    // Judd(1998), page 233,   Algorithm 6.3 along with equations 6.11.4 and
    // 6.11.5 provide this whole section.
    const auto delta = (y[1] - y[0]) / (x[1] - x[0]);
    const auto Condition =
        ((gradients[0] - delta) * (gradients[1] - delta) >= 0);
    const auto Condition2 =
        std::abs(gradients[1] - delta) < std::abs(gradients[0] - delta);
    if (Condition) {
      tsi = sum_vector(x) / 2;
    } else if (Condition2) {
      tsi = (x[0] + (x[1] - x[0]) * (gradients[1] - delta) /
                        (gradients[1] - gradients[0]));
    } else {
      tsi = (x[1] + (x[1] - x[0]) * (gradients[0] - delta) /
                        (gradients[1] - gradients[0]));
    }
  }

  // #Judd(1998), page 232, 3rd last equation of page.
  const auto alpha = tsi - x[0];
  const auto beta = x[1] - tsi;
  // #Judd(1998), page 232, 4th last equation of page.
  const auto sbar =
      (2 * (y[1] - y[0]) - (alpha * gradients[0] + beta * gradients[1])) /
      (x[1] - x[0]);

  // #Judd(1998), page 232, 3rd equation of page.(C1, B1, A1)
  std::vector<double> Coeffs1 = {(sbar - gradients[0]) / (2 * alpha),
                                 gradients[0], y[0]};
  std::vector<double> Coeffs2;
  if (beta == 0) {
    Coeffs2 = Coeffs1;
  } else {
    // #Judd(1998), page 232, 4th equation of page.(C2, B2, A2)
    double vals = 0;
    for (int i = 0; i < 3; i++) {
      vals += Coeffs1[i] * std::pow(alpha, 2 - i);
    }
    Coeffs2 = {(gradients[1] - sbar) / (2 * beta), sbar, vals};
  }

  // returning the answer
  const auto Machine4Epsilon = 4 * DBL_EPSILON;
  std::vector<double> intervals;
  std::vector<std::vector<double>> result;
  if (tsi < x[0] + Machine4Epsilon) {
    intervals = {x[0]};
    result.push_back(Coeffs2);
  } else if (tsi > x[1] - Machine4Epsilon) {
    intervals = {x[0]};
    result.push_back(Coeffs1);
  } else {
    intervals = {x[0], tsi};
    result.push_back(Coeffs1);
    result.push_back(Coeffs2);
  }
  return {intervals, result};
}

/*
 Calls `SchumakerIndInterval` many times to get full set of spline intervals
and coefficients. Then calls extrapolation for out of sample behaviour
### Inputs
 * `gradients` - A vector of gradients at each point
 * `x` - A vector of `x` coordinates
 * `y` - A vector of `y` coordinates
 * `extrapolation` - A string in (Curve, Linear or Constant) that gives
behaviour outside of interpolation range.
### Returns
 * A vector of interval starts
 * A vector of interval ends
 * A matrix of all coefficients
*/
std::tuple<std::vector<double>, std::vector<std::vector<double>>>
getCoefficientMatrix(const std::vector<double> x, const std::vector<double> y,
                     const std::vector<double> gradients) {
  const auto n = x.size();
  std::vector<double> allstarts;
  std::vector<std::vector<double>> fullMatrix;
  for (int intrval = 1; intrval < n; intrval++) {
    const std::vector<double> xs = {x[intrval - 1], x[intrval]};
    const std::vector<double> ys = {y[intrval - 1], y[intrval]};
    const std::vector<double> grads = {gradients[intrval - 1],
                                       gradients[intrval]};
    auto [starts, intMatrix] = schumakerIndInterval(grads, ys, xs);
    allstarts = append_vectors(allstarts, starts);
    fullMatrix = append_vectors(fullMatrix, intMatrix);
  }
  return {allstarts, fullMatrix};
}

/*
  imputeGradients(x::Vector{T}, y::Vector{T})
Imputes gradients based on a vector of x and y coordinates.
  */
std::vector<double> imputeGradients(const std::vector<double> x,
                                    const std::vector<double> y) {
  const auto n = x.size();
  std::vector<double> L(n);
  std::vector<double> d(n);
  for (int i = 1; i < n + 1; i++) {
    // Judd (1998), page 233, second last equation
    L[i - 1] =
        std::sqrt(std::pow(x[i] - x[i - 1], 2) + std::pow(y[i] - y[i - 1], 2));
    // Judd (1998), page 233, last equation
    d[i - 1] = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
  }

  std::vector<double> sb(n);
  for (int i = 1; i < n + 1; i++) {
    // Judd (1998), page 234, Eqn 6.11.6
    const auto cond = d[i - 1] * d[i] > 0;
    const auto val = (L[i - 1] * d[i - 1] + L[i] * d[i]) / (L[i - 1] + L[i]);
    sb[i - 1] = cond * val;
  }

  // Judd (1998), page 234, Second Equation line plus 6.11.6 gives this array
  // of slopes.
  std::vector<double> grads(n);
  grads[0] = ((-sb[0] + 3 * d[0]) / 2);
  for (auto i = 1; i < n - 1; i++) {
    grads[i] = sb[i - 1];
  }
  grads[n - 1] = ((3 * d[n - 2] - sb[n - 3]) / 2);
  return grads;
}

Spline::Spline(const std::vector<double> x,
               const std::vector<double> y) { // Constructor with dataset
  if (x.size() != y.size()) {
    throw std::invalid_argument("x and y vectors must have the same size.");
  }
  const auto &gradients = imputeGradients(x, y);
  const auto [intervals, matrix] = getCoefficientMatrix(x, y, gradients);
  this->Intervals_ = intervals;
  this->Coefficients_ = matrix;
}

double Spline::spline(
    double inp) { // This gives you the value of the spline at the input.
  const auto &ind = search_sorted_last_less_than_x(this->Intervals_, inp);
  const auto &Coeffs = this->Coefficients_[ind];
  const auto &rescaled_inp = inp - this->Intervals_[ind];

  return Coeffs[0] * std::pow(rescaled_inp, 2) +
         Coeffs[1] * std::pow(rescaled_inp, 1) + Coeffs[2];
}

double Spline::derivative(
    double inp) { // This gives you the value of the derivative at the input
  const auto &ind = search_sorted_last_less_than_x(this->Intervals_, inp);
  const auto &Coeffs = this->Coefficients_[ind];
  const auto &rescaled_inp = inp - this->Intervals_[ind];
  return 2 * Coeffs[0] * std::pow(rescaled_inp, 1) + Coeffs[1];
}

double Spline::second_derivative(double inp) { // This gives you the value of the
                                        // second derivative at the input
  const auto &ind = search_sorted_last_less_than_x(this->Intervals_, inp);
  const auto &Coeffs = this->Coefficients_[ind];
  return 2 * Coeffs[0];
}

std::vector<double> Spline::spline(
    std::vector<double> inp) { // This gives you the value of the spline at the input.
  auto result = std::vector<double>(inp.size());
  for (auto i = 0; i < inp.size(); i++) {
    result[i] = spline(inp[i]);
  }
  return result;
}

std::vector<double> Spline::derivative(
    std::vector<double> inp) { // This gives you the value of the derivative at the input
  auto result = std::vector<double>(inp.size());
  for (auto i = 0; i < inp.size(); i++) {
    result[i] = derivative(inp[i]);
  }
  return result;
}

std::vector<double> Spline::second_derivative(std::vector<double> inp) { // This gives you the value of the
                                        // second derivative at the input
  auto result = std::vector<double>(inp.size());
  for (auto i = 0; i < inp.size(); i++) {
    result[i] = second_derivative(inp[i]);
  }
  return result;
}
