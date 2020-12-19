//
// Created by Avinash Madavan on 2019-06-28.
//

#ifndef STUKA_LP_PRECONDITIONED_LP_H
#define STUKA_LP_PRECONDITIONED_LP_H

#include "base_lp.h"
#include "preconditioner/base_preconditioner.h"

namespace stuka { namespace LP {

template<class LP, class Preconditioner>
class PreconditionedLinearProgram : public BaseLinearProgram {
 public:
  explicit PreconditionedLinearProgram() : lp(std::make_shared<LP>()) {
    std::shared_ptr<BaseLinearProgram> lp_base = lp;
    next_ = std::make_shared<Preconditioner>(lp_base);
  }

  ~PreconditionedLinearProgram() override = default;

  void initialize(const LinearProgram &prog) override { next_->initialize(prog); }

  void setObjective(const std::unique_ptr<Eigen::VectorXd> &c) override { next_->setObjective(c); }

  void setRHS(const std::unique_ptr<Eigen::VectorXd> &b_ub, const std::unique_ptr<Eigen::VectorXd> &b_eq) override {
    next_->setRHS(b_ub, b_eq);
  }

  void setBounds(const std::unique_ptr<Eigen::VectorXd> &lb, const std::unique_ptr<Eigen::VectorXd> &ub) override {
    next_->setBounds(lb, ub);
  }

  void addVar(double c, const std::unique_ptr<Eigen::VectorXd> &a_ub, const std::unique_ptr<Eigen::VectorXd> &a_eq, double lb,
              double ub) override { next_->addVar(c, a_ub, a_eq, lb, ub); }

  void addVars(const std::unique_ptr<Eigen::VectorXd> &c, const std::unique_ptr<Eigen::SparseMatrix<double>> &A_ub,
               const std::unique_ptr<Eigen::SparseMatrix<double>> &A_eq, const std::unique_ptr<Eigen::VectorXd> &lb,
               const std::unique_ptr<Eigen::VectorXd> &ub) override { next_->addVars(c, A_ub, A_eq, lb, ub); }

  void removeVar(size_t var) override { next_->removeVar(var); }

  void removeVars(size_t index, size_t n_remove) override { next_->removeVars(index, n_remove); }

  void removeBackVars(size_t n_remove) override { next_->removeBackVars(n_remove); }

  void addConstr_ub(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) override { next_->addConstr_ub(a, b); }

  void addConstrs_ub(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                     const std::unique_ptr<Eigen::VectorXd> &b) override { next_->addConstrs_ub(A, b); }

  void removeConstr_ub(size_t index) override { next_->removeConstr_ub(index); }

  void removeConstrs_ub(size_t index, size_t n_remove) override { next_->removeConstrs_ub(index, n_remove); }

  void addConstr_eq(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) override { next_->addConstr_eq(a, b); }

  void addConstrs_eq(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                     const std::unique_ptr<Eigen::VectorXd> &b) override { next_->addConstrs_eq(A, b); }

  void removeConstr_eq(size_t index) override { next_->removeConstr_eq(index); }

  void removeConstrs_eq(size_t index, size_t n_remove) override { next_->removeConstrs_eq(index, n_remove); }

  Eigen::VectorXd convertState(const Eigen::VectorXd &x) override { return next_->convertState(x); }

  Eigen::VectorXd revertState(const Eigen::VectorXd &x) override { return next_->revertState(x); }

  const std::shared_ptr<LP> lp;

 private:
  std::shared_ptr<BaseLinearProgram> next_;
};
}}

#endif //STUKA_LP_PRECONDITIONED_LP_H
