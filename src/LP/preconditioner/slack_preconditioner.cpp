//
// Created by Avinash Madavan on 2019-06-26.
//

#include <stuka/LP/preconditioner/slack_preconditioner.h>

stuka::LP::SlackPreconditioner::SlackPreconditioner(const std::shared_ptr<stuka::LP::BaseLinearProgram> &next)
    : BasePreconditioner(next) {}

void stuka::LP::SlackPreconditioner::initialize(const stuka::LP::LinearProgram &prog) {
  n_dim_ = prog.c->size();
  n_ub_ = (prog.b_ub) ? prog.b_ub->size(): 0;
  n_eq_ = (prog.b_eq) ? prog.b_eq->size(): 0;

  A_ub_ = prog.A_ub;
  b_ub_ = prog.b_ub;
  A_eq_ = prog.A_eq;
  b_eq_ = prog.b_eq;
//  A_slack_ = std::make_shared<Eigen::

  if (n_ub_ == 0) {
    next()->initialize(prog);
  }
  else {
    LinearProgram lp_next;
    lp_next.c = std::make_shared<Eigen::VectorXd>(n_dim_ + n_ub_);
    lp_next.c->head(n_dim_) = *prog.c;

    lp_next.A_eq = std::make_shared<Eigen::SparseMatrix<double>>(n_eq_ + n_ub_, n_dim_ + n_ub_);
    lp_next.A_eq->reserve(n_ub_ + ((n_eq_ > 0) ? prog.A_eq->nonZeros() : 0) + ((n_ub_ > 0) ? prog.A_ub->nonZeros() : 0));
    for (size_t i = 0; i < n_dim_; ++i) {
      lp_next.A_eq->startVec(i);
      if (n_eq_ > 0)
        for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_eq, i); it; ++it)
          lp_next.A_eq->insertBack(it.row(), i) = it.value();
      for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_ub, i); it; ++it)
        lp_next.A_eq->insertBack(n_eq_ + it.row(), i) = it.value();
    }

    for (size_t i = 0; i < n_ub_; ++i) {
      lp_next.A_eq->startVec(n_dim_ + i);
      lp_next.A_eq->insertBack(n_eq_ + i, n_dim_ + i) = 1.;
    }

    lp_next.A_eq->finalize();

    lp_next.b_eq = std::make_shared<Eigen::VectorXd>(n_eq_ + n_ub_);
    if (n_eq_ > 0) lp_next.b_eq->head(n_eq_) = *prog.b_eq;
    lp_next.b_eq->tail(n_ub_) = *prog.b_ub;

    lp_next.lb = std::make_shared<Eigen::VectorXd>(n_dim_ + n_ub_);
    lp_next.lb->setZero();
    if (prog.lb) lp_next.lb->head(n_dim_) = *prog.lb;

    if (prog.ub) {
      lp_next.ub = std::make_shared<Eigen::VectorXd>(n_dim_ + n_ub_);
      lp_next.ub->setConstant(INF);
      lp_next.ub->head(n_dim_) = *prog.ub;
    }

    next()->initialize(lp_next);
  }
}

void stuka::LP::SlackPreconditioner::setObjective(const std::shared_ptr<Eigen::VectorXd> &c) {
  std::shared_ptr<Eigen::VectorXd> c_next = std::make_shared<Eigen::VectorXd>(*c);
  c_next->conservativeResize(n_dim_ + n_ub_);
  c_next->tail(n_ub_).setZero();
  next()->setObjective(c_next);
}

void stuka::LP::SlackPreconditioner::setRHS(const std::shared_ptr<Eigen::VectorXd> &b_ub,
                                            const std::shared_ptr<Eigen::VectorXd> &b_eq) {
  if (b_ub) b_ub_ = b_ub;
  if (b_eq) b_eq_ = b_eq;

  std::shared_ptr<Eigen::VectorXd> b_eq_next = std::make_shared<Eigen::VectorXd>(n_eq_ + n_ub_);
  b_eq_next->head(n_eq_) = *b_eq_;
  b_eq_next->tail(n_ub_) = *b_ub_;

  next()->setRHS(nullptr, b_eq_next);
}

void stuka::LP::SlackPreconditioner::setBounds(const std::shared_ptr<Eigen::VectorXd> &lb,
                                               const std::shared_ptr<Eigen::VectorXd> &ub) {
  next()->setBounds(lb, ub);
}

void stuka::LP::SlackPreconditioner::addVar(double c, std::shared_ptr<Eigen::VectorXd> a_ub,
                                            std::shared_ptr<Eigen::VectorXd> a_eq, double lb, double ub) {

}

void stuka::LP::SlackPreconditioner::addVars(std::shared_ptr<Eigen::VectorXd> c,
                                             std::shared_ptr<Eigen::SparseMatrix<double>> A_ub,
                                             std::shared_ptr<Eigen::SparseMatrix<double>> A_eq,
                                             std::shared_ptr<Eigen::VectorXd> lb, std::shared_ptr<Eigen::VectorXd> ub) {

}

void stuka::LP::SlackPreconditioner::removeVar(size_t var) {

}

void stuka::LP::SlackPreconditioner::removeVars(size_t index, size_t n_remove) {

}

void stuka::LP::SlackPreconditioner::removeBackVars(size_t n_remove) {

}

void stuka::LP::SlackPreconditioner::addConstr_ub(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) {

}

void stuka::LP::SlackPreconditioner::addConstrs_ub(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                                                   const std::shared_ptr<Eigen::VectorXd> &b) {

}

void stuka::LP::SlackPreconditioner::removeConstr_ub(size_t index) {

}

void stuka::LP::SlackPreconditioner::removeConstrs_ub(size_t index, size_t n_remove) {

}

void stuka::LP::SlackPreconditioner::addConstr_eq(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) {

}

void stuka::LP::SlackPreconditioner::addConstrs_eq(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                                                   const std::shared_ptr<Eigen::VectorXd> &b) {

}

void stuka::LP::SlackPreconditioner::removeConstr_eq(size_t index) {

}

void stuka::LP::SlackPreconditioner::removeConstrs_eq(size_t index, size_t n_remove) {

}

Eigen::VectorXd stuka::LP::SlackPreconditioner::convertState(const Eigen::VectorXd &x) {
  Eigen::VectorXd x_(n_dim_ + n_ub_);
  x_.setZero();
  x_.head(n_dim_) = x;
  x_.tail(n_ub_) = *b_ub_ - *A_ub_ * x;
  return next()->convertState(x_);
}

Eigen::VectorXd stuka::LP::SlackPreconditioner::revertState(const Eigen::VectorXd &x) {
  return next()->revertState(x).head(n_dim_);
}
