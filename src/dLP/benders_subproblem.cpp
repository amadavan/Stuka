//
// Created by Avinash Madavan on 1/15/19.
//

#include <stuka/dLP/benders_subproblem.h>

stuka::dLP::BendersSubproblem::BendersSubproblem(const stuka::dLP::Subproblem sub, const Options opts) : sub_(sub) {
  current_value_ = INF;

  LP::LinearProgram lp;
  lp.c = sub.c;
  lp.A_ub = sub.A_ub;
  lp.b_ub = sub.b_ub;
  lp.A_eq = sub.A_eq;
  lp.b_eq = sub.b_eq;
  lp.lb = sub.lb;
  lp.ub = sub.ub;

  solver_ = util::createSolver(lp, opts);
}

stuka::dLP::BendersSubproblem::BendersSubproblem(stuka::dLP::BendersSubproblem &&sub) noexcept :
    sub_(sub.sub_), solver_(std::move(sub.solver_)) {
  current_value_ = INF;
}

stuka::dLP::BendersCut stuka::dLP::BendersSubproblem::getBendersCut(const Eigen::VectorXd &x) {

  solver_->getLP().setRHS((sub_.b_ub) ? std::make_shared<Eigen::VectorXd>(*sub_.b_ub - *sub_.C_ub * x) : nullptr,
                          (sub_.b_eq) ? std::make_shared<Eigen::VectorXd>(*sub_.b_eq - *sub_.C_eq * x) : nullptr);

  // Solve the subproblem at the given point
  OptimizeState res;
  try {
    res = solver_->solve();
  } catch (std::exception &e) {
    throw std::runtime_error("getBendersCut: unable to solve subproblem");
  }
  if (res.status != 2) {
    std::cout << "Failed to solve subproblem: " << res.status << std::endl;
    throw std::runtime_error("getBendersCut: unable to solve subproblem");
  }

  Eigen::VectorXd aux = Eigen::VectorXd::Zero(x.size());
  if (sub_.b_ub)
    aux += sub_.C_ub->transpose() * res.dual_ub;
  if (sub_.b_eq)
    aux += sub_.C_eq->transpose() * res.dual_eq;

  current_value_ = res.fun;

  return {aux, -res.fun + aux.dot(x)};
}

const double &stuka::dLP::BendersSubproblem::getValue() {
  return current_value_;
}
