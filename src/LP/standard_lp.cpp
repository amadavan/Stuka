//
// Created by Avinash Madavan on 2019-07-03.
//

#include <stuka/LP/standard_lp.h>

void stuka::LP::StandardLinearProgram::initialize(const stuka::LP::LinearProgram &prog) {
  c = util::DenseOps::unique_copy(prog.c);
  A_ub = util::SparseOps::unique_copy(prog.A_ub);
  b_ub = util::DenseOps::unique_copy(prog.b_ub);
  A_eq = util::SparseOps::unique_copy(prog.A_eq);
  b_eq = util::DenseOps::unique_copy(prog.b_eq);
  lb = util::DenseOps::unique_copy(prog.lb);
  ub = util::DenseOps::unique_copy(prog.ub);

  n_dim_ = this->c->size();
  n_ub_ = (this->b_ub) ? this->b_ub->size() : 0;
  n_eq_ = (this->b_eq) ? this->b_eq->size() : 0;
}

void stuka::LP::StandardLinearProgram::setObjective(const std::unique_ptr<Eigen::VectorXd> &c) {
  this->c = util::DenseOps::unique_copy(c);
}

void stuka::LP::StandardLinearProgram::setRHS(const std::unique_ptr<Eigen::VectorXd> &b_ub,
                                              const std::unique_ptr<Eigen::VectorXd> &b_eq) {
  this->b_ub = util::DenseOps::unique_copy(b_ub);
  this->b_eq = util::DenseOps::unique_copy(b_eq);
}

void stuka::LP::StandardLinearProgram::setBounds(const std::unique_ptr<Eigen::VectorXd> &lb,
                                                 const std::unique_ptr<Eigen::VectorXd> &ub) {
  this->lb = util::DenseOps::unique_copy(lb);
  this->ub = util::DenseOps::unique_copy(ub);
}

void stuka::LP::StandardLinearProgram::addVar(double c, const std::unique_ptr<Eigen::VectorXd> &a_ub,
                                              const std::unique_ptr<Eigen::VectorXd> &a_eq, double lb, double ub) {
  this->c->conservativeResize(n_dim_ + 1);
  this->c->coeffRef(n_dim_ + 1) = c;

  // TODO: other stuff
}

void stuka::LP::StandardLinearProgram::addVars(const std::unique_ptr<Eigen::VectorXd> &c,
                                               const std::unique_ptr<Eigen::SparseMatrix<double>> &A_ub,
                                               const std::unique_ptr<Eigen::SparseMatrix<double>> &A_eq,
                                               const std::unique_ptr<Eigen::VectorXd> &lb,
                                               const std::unique_ptr<Eigen::VectorXd> &ub) {

}

void stuka::LP::StandardLinearProgram::removeVar(size_t var) {

}

void stuka::LP::StandardLinearProgram::removeVars(size_t index, size_t n_remove) {

}

void stuka::LP::StandardLinearProgram::removeBackVars(size_t n_remove) {

}

void stuka::LP::StandardLinearProgram::addConstr_ub(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) {

}

void stuka::LP::StandardLinearProgram::addConstrs_ub(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                                                     const std::unique_ptr<Eigen::VectorXd> &b) {

}

void stuka::LP::StandardLinearProgram::removeConstr_ub(size_t index) {

}

void stuka::LP::StandardLinearProgram::removeConstrs_ub(size_t index, size_t n_remove) {

}

void stuka::LP::StandardLinearProgram::addConstr_eq(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) {

}

void stuka::LP::StandardLinearProgram::addConstrs_eq(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                                                     const std::unique_ptr<Eigen::VectorXd> &b) {

}

void stuka::LP::StandardLinearProgram::removeConstr_eq(size_t index) {

}

void stuka::LP::StandardLinearProgram::removeConstrs_eq(size_t index, size_t n_remove) {

}

Eigen::VectorXd stuka::LP::StandardLinearProgram::convertState(const Eigen::VectorXd &x) {
  return x;
}

Eigen::VectorXd stuka::LP::StandardLinearProgram::revertState(const Eigen::VectorXd &x) {
  return x;
}
