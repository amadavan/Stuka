//
// Created by Avinash Madavan on 2019-01-29.
//

#ifndef STUKA_UTIL_FUNCTIONS_H
#define STUKA_UTIL_FUNCTIONS_H

#include <memory>

#include "../LP/lp.h"
#include "../LP/base_solver.h"
#include "../QP/qp.h"
#include "../QP/base_solver.h"
#include "../dLP/decomposed_lp.h"
#include "../dLP/base_solver.h"
#include "../stochastic/program.h"

#include "../optimize_state.h"
#include "../options.h"

#include "timer.h"
#include "solver_factory.h"

namespace stuka { namespace util {
const OptimizeState linprog(const LP::LinearProgram &lp, const Options &opts = Options());

const OptimizeState quadprog(const QP::QuadraticProgram &qp, const Options &opts = Options());

const OptimizeState linprog(const dLP::DecomposedLinearProgram &lp, const Options &opts = Options());

const OptimizeState stochastic_cvar(const stochastic::Program &prog, const Options &opts = Options());
}}

#endif //STUKA_UTIL_FUNCTIONS_H
