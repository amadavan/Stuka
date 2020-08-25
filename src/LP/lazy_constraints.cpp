//
// Created by Avinash on 8/19/2020.
//

#include <stuka/LP/lazy_constraints.h>

stuka::LP::LazyConstraints::LazyConstraints(const stuka::LP::LinearProgram &lp, const stuka::Options &opts)
    : BaseLPSolver(lp, opts), lp_() {

  stuka::Options _opts(opts);
  _opts.lazy = false;

  LinearProgram lp2(lp);
  lp2.A_ub = nullptr;
  lp2.b_ub = nullptr;
  solver_ = util::createSolver(lp2, _opts);

  lp_.setLP(&solver_->getLP());
  lp_.initialize(lp);

  nit = 0;
}

void stuka::LP::LazyConstraints::iterate() {
  nit++;

  try {
    state_ = solver_->solve();
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    throw std::runtime_error("Lazy constraint generation had a problem with gurobi...");
  }

  // TODO: Determine a better way to evaluate constraint generation for unbounded solutions.
  if (state_.status != 2) lp_.addRandomConstraint(); // If unbounded, then constraints need to be added.

  lp_.addViolatedConstraints(state_.x);
}

bool stuka::LP::LazyConstraints::terminate() {
  // TODO: This is a very brittle termination criterion. Should verify that all constraints have been added instead.
  return (state_.x.size() == 0 || !lp_.isViolated(state_.x)) && state_.status == 2;
}

const stuka::OptimizeState stuka::LP::LazyConstraints::getState() {
  OptimizeState state(state_);
  if (state.status == 2) state.dual_ub = lp_.getDualUB(state.dual_ub);
  state.nit_con_gen = nit;
  return state;
}

stuka::LP::BaseLinearProgram &stuka::LP::LazyConstraints::getLP() {
  return lp_;
}
