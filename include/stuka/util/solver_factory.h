//
// Created by Avinash Madavan on 1/16/19.
//

#ifndef STUKA_UTIL_SOLVER_FACTORY_H
#define STUKA_UTIL_SOLVER_FACTORY_H

#include <memory>

#include "../LP/lp.h"
#include "../LP/base_solver.h"
#include "../LP/lazy_constraints.h"
#include "../LP/mehrotra_pc.h"
#include "../QP/qp.h"
#include "../QP/base_solver.h"
#include "../dLP/decomposed_lp.h"
#include "../dLP/base_solver.h"
#include "../stochastic/base_solver.h"
#include "../stochastic/stochastic_pd.h"
#include "../stochastic/stochastic_pd2.h"
#include "../stochastic/stochastic_pd_jacobi.h"

#ifdef ENABLE_GUROBI
#include "../gurobi/LP/gurobi.h"
#include "../gurobi/QP/gurobi.h"
#endif

namespace stuka { namespace util {
std::unique_ptr<LP::BaseLPSolver> createSolver(const LP::LinearProgram &lp, const Options &opts = Options());

std::unique_ptr<QP::BaseQPSolver> createSolver(const QP::QuadraticProgram &qp, const Options &opts = Options());

std::unique_ptr<dLP::BaseDLPSolver>
createSolver(const dLP::DecomposedLinearProgram &dlp, const Options &opts = Options());

std::unique_ptr<stochastic::BaseStochasticSolver>
createSolver(const stochastic::Program &prog, const Options &opts = Options());
}}
#endif //STUKA_UTIL_SOLVER_FACTORY_H
