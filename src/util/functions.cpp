//
// Created by Avinash Madavan on 2019-01-29.
//

#include <stuka/util/functions.h>

const stuka::OptimizeState stuka::util::linprog(const stuka::LP::LinearProgram &lp, const stuka::Options &opts) {
  Timer timer;

  timer.start();
  std::unique_ptr<LP::BaseLPSolver> solver = createSolver(lp, opts);
  OptimizeState res = solver->solve();
  res.runtime = timer.elapsed().count();

  return res;
}

const stuka::OptimizeState stuka::util::quadprog(const stuka::QP::QuadraticProgram &qp, const stuka::Options &opts) {
  Timer timer;

  timer.start();
  std::unique_ptr<QP::BaseQPSolver> solver = createSolver(qp, opts);
  OptimizeState res = solver->solve();
  res.runtime = timer.elapsed().count();

  return res;
}

const stuka::OptimizeState
stuka::util::linprog(const stuka::dLP::DecomposedLinearProgram &dlp, const stuka::Options &opts) {
  Timer timer;

  timer.start();
  std::unique_ptr<dLP::BaseDLPSolver> solver = createSolver(dlp, opts);
  OptimizeState res = solver->solve();
  res.runtime = timer.elapsed().count();

  return res;
}

const stuka::OptimizeState
stuka::util::stochastic_cvar(const stuka::stochastic::Program &prog, const stuka::Options &opts) {
  Timer timer;

  timer.start();
  std::unique_ptr<stochastic::BaseStochasticSolver> solver = createSolver(prog, opts);
  OptimizeState res = solver->solve();
  res.runtime = timer.elapsed().count();

  return res;
}
