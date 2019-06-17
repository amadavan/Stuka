//
// Created by Avinash Madavan on 1/16/19.
//

#include <stuka/util/solver_factory.h>
#include <stuka/stuka.h>

std::unique_ptr<stuka::LP::BaseLPSolver>
stuka::util::createSolver(const LP::LinearProgram &lp, const stuka::Options &opts) {
  std::unique_ptr<LP::BaseLPSolver> solver;

  switch (opts.lp_solver) {
#ifdef ENABLE_GUROBI
    case GUROBI:
      solver = std::make_unique<LP::GurobiSolver>(lp, opts); break;
#endif
    case MPC:
      solver = std::make_unique<LP::MehrotraPC>(lp, opts);
      break;
    default:
      throw std::runtime_error("createSolver: invalid solver selected");
  }

  return solver;
}

std::unique_ptr<stuka::QP::BaseQPSolver>
stuka::util::createSolver(const stuka::QP::QuadraticProgram &qp, const stuka::Options &opts) {
  std::unique_ptr<QP::BaseQPSolver> solver;

  switch (opts.qp_solver) {
#ifdef ENABLE_GUROBI
    case GUROBI:
      solver = std::make_unique<QP::GurobiSolver>(qp, opts); break;
#endif
    default:
      throw std::runtime_error("createSolver: invalid solver selected");
  }

  return solver;
}

std::unique_ptr<stuka::dLP::BaseDLPSolver>
stuka::util::createSolver(const stuka::dLP::DecomposedLinearProgram &dlp, const stuka::Options &opts) {
  std::unique_ptr<dLP::BaseDLPSolver> solver;

  switch (opts.dlp_solver) {
    case BENDER:
      solver = std::make_unique<dLP::BendersDecomposition>(dlp, opts); break;
    case CRE:
      solver = std::make_unique<dLP::CRE>(dlp, opts); break;
    default:
      throw std::runtime_error("createSolver: invalid solver selected");
  }

  return solver;
}
