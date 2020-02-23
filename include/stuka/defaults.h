//
// Created by Avinash Madavan on 2/12/20.
//

#ifndef STUKA_DEFAULTS_H_
#define STUKA_DEFAULTS_H_

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

// Define default cre step size
#ifndef DEFAULT_CRE_STEP
#define DEFAULT_CRE_STEP 1e-4
#endif

#endif //STUKA_DEFAULTS_H_
