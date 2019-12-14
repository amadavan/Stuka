//
// Created by Avinash Madavan on 1/7/19.
//

#ifndef STUKA_CONSTANTS_H
#define STUKA_CONSTANTS_H

#include <stddef.h>

#include <cmath>

#include "util/common_expressions.h"

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
};

enum ConstraintReductionMethods {
  BOUNDS
};

enum Status {

};

}

#endif //STUKA_CONSTANTS_H
