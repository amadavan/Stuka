//
// Created by Avinash Madavan on 1/7/19.
//

#ifndef STUKA_CONSTANTS_H
#define STUKA_CONSTANTS_H

#include <stddef.h>

#include <cmath>

#include "util/common_expressions.h"

// Define Gurobi options
#ifdef ENABLE_GUROBI
#ifndef GUROBI_TOLERANCE
#define GUROBI_TOLERANCE 1e-9
#endif
#endif

#include <limits>

namespace stuka {

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
  PDSS_JACOBI,
  NAIVE_BENDERS,
  NAIVE_CRE
};

enum ConstraintReductionMethods {
  BOUNDS
};

enum Status {

};

}

#endif //STUKA_CONSTANTS_H
