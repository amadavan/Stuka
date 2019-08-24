//
// Created by Avinash Madavan on 1/7/19.
//

#include <stuka/QP/gurobi.h>

stuka::QP::GurobiSolver::GurobiSolver(const stuka::QP::QuadraticProgram &qp, const stuka::Options &opts) : BaseQPSolver(
    qp, opts), prog_(qp) {}

stuka::QP::BaseQuadraticProgram &stuka::QP::GurobiSolver::getQP() {
  return prog_;
}

void stuka::QP::GurobiSolver::iterate() {
  prog_.model_.optimize();
}

bool stuka::QP::GurobiSolver::terminate() {
  return true;
}

const stuka::OptimizeState stuka::QP::GurobiSolver::getState() {
  OptimizeState state;
  state.x = Eigen::VectorXd(prog_.n_dim_);

  state.status = prog_.model_.get(GRB_IntAttr_Status);

  if (state.status == 2) {
    state.fun = prog_.model_.get(GRB_DoubleAttr_ObjVal);

    state.x = Eigen::Map<Eigen::VectorXd>(prog_.model_.get(GRB_DoubleAttr_X, prog_.vars_, (int) prog_.n_dim_),
                                          prog_.n_dim_);

    if (prog_.n_con_ub_ > 0) {
      state.dual_ub = Eigen::Map<Eigen::VectorXd>(
          prog_.model_.get(GRB_DoubleAttr_Pi, prog_.ubconstr_, (int) prog_.n_con_ub_), prog_.n_con_ub_);
      state.dual_ub *= -1;
    }

    if (prog_.n_con_eq_ > 0) {
      state.dual_eq = Eigen::Map<Eigen::VectorXd>(
          prog_.model_.get(GRB_DoubleAttr_Pi, prog_.eqconstr_, (int) prog_.n_con_eq_), prog_.n_con_eq_);
      state.dual_eq *= -1;
    }

    Eigen::VectorXd rc = Eigen::Map<Eigen::VectorXd>(prog_.model_.get(GRB_DoubleAttr_RC, prog_.model_.getVars(), (int) prog_.n_dim_), prog_.n_dim_);
    state.dual_x_lb = Eigen::VectorXd(prog_.n_dim_);
    state.dual_x_ub = Eigen::VectorXd(prog_.n_dim_);

    state.dual_x_lb.setZero();
    state.dual_x_ub.setZero();

    for (size_t i = 0; i < prog_.n_dim_; ++i) {
      if (rc[i] < -1e-9)
        state.dual_x_ub.coeffRef(i) = -rc.coeff(i);
      else if (rc[i] > 1e-9)
        state.dual_x_lb.coeffRef(i) = rc.coeff(i);
    }
  } else {
    state.fun = INF;
  }

  return state;
}
