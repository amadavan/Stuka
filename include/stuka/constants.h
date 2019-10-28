//
// Created by Avinash Madavan on 1/7/19.
//

#ifndef STUKA_CONSTANTS_H
#define STUKA_CONSTANTS_H

#define _USE_MATH_DEFINES

#include <stddef.h>

#include <cmath>

// Define default LP solver
#ifndef DEFAULT_LP_SOLVER
#ifdef ENABLE_GUROBI
#define DEFAULT_LP_SOLVER GUROBI
#else
#define DEFAULT_LP_SOLVER MPC
#endif
#endif

// Define default QP solver
#ifndef DEFAULT_QP_SOLVER
#ifdef ENABLE_GUROBI
#define DEFAULT_QP_SOLVER GUROBI
#else
#define DEFAULT_QP_SOLVER MPC
#endif
#endif

// Define default dLP solver
#ifndef DEFAULT_DLP_SOLVER
#define DEFAULT_DLP_SOLVER BENDER
#endif

// Define default stochastic solver
#ifndef DEFAULT_STOCHASTIC_SOLVER
#define DEFAULT_STOCHASTIC_SOLVER PDSS
#endif

#include <limits>

#define EXP_MAX_RECURSION_DEPTH 64
#define PI_MAX_RECURSION_DEPTH 64
#define SQRT_MAX_RECURSION_DEPTH 500

namespace stuka {
template<typename T, typename P, std::enable_if_t<std::is_integral<P>::value, int> = 0>
constexpr T pow(T x, P y) {
  return y == 0 ? 1.0 : x * pow(x, y - 1);
}

// Estimate the exponential function using a Taylor series expansion
template<typename T>
constexpr T exp_recursion(T x, T factorial_inverse, const size_t depth) noexcept {
  if (depth < EXP_MAX_RECURSION_DEPTH) {
    factorial_inverse *= (depth == 0) ? 1 : x / depth;
    return factorial_inverse + exp_recursion<T>(x, factorial_inverse, depth + 1);
  } else {
    return 0;
  }
}

template<typename T>
constexpr T exp(T x) {
  return exp_recursion<T>(x, 1, 0);
}

// Estimating pi is more difficult, so we just use the actual value
template<typename T>
constexpr T pi_est() {
  return T(3.14159265358979323846264338327950288419716939937510582097494459230781640628L);
}

// Estimate pi with Chudnovsky algorithm
template<typename T>
constexpr T chudnovsky_recursion(T L, T X, T K, T M, size_t depth) {
  if (depth < PI_MAX_RECURSION_DEPTH) {
    T current_value = M * L / X;
    L += 545140134;
    X *= -262537412640768000;
    M *= (pow(K, 3) - 16 * K) / pow(T(depth) + 1, 3);
    K += 12;

    return current_value + chudnovsky_recursion(L, X, K, M, depth + 1);
  } else {
    return 0;
  }
}

template<typename T>
constexpr T pi_chudnovsky() {
  T sqrt_10005 = 100.0249968757810059447921878763577780015950243686963146571355115696509678538643042923111879484999733L;
  return 426880 * sqrt_10005 / chudnovsky_recursion<T>(13591409, 1, 6, 1, 0);
}

constexpr double INF = std::numeric_limits<double>::infinity();
constexpr long double e = exp<long double>(1.);
constexpr long double pi = pi_chudnovsky<long double>();

enum Solver {
  GUROBI,
  BENDER,
  CRE,
  MPC,
  PDSS,
  PDSS2,
};

enum ConstraintReductionMethods {
  BOUNDS
};

enum Status {

};

}

#endif //STUKA_CONSTANTS_H
