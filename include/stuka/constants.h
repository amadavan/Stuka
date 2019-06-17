//
// Created by Avinash Madavan on 1/7/19.
//

#ifndef STUKA_CONSTANTS_H
#define STUKA_CONSTANTS_H

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

#include <limits>

namespace stuka {

  constexpr double INF = std::numeric_limits<double>::infinity();

  enum Solver {
#ifdef ENABLE_GUROBI
    GUROBI,
#endif
    BENDER,
    CRE,
    MPC
  };

}

#endif //STUKA_CONSTANTS_H
