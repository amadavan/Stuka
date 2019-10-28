//
// Created by Avinash Madavan on 2019-07-10.
//

#ifndef STUKA_LP_COMPOSITE_PRECONDITIONER_H
#define STUKA_LP_COMPOSITE_PRECONDITIONER_H

#include "../base_lp.h"
#include "base_preconditioner.h"

namespace stuka { namespace LP {

template<class ...Preconditioners>
class CompositePreconditioner : public BasePreconditioner {
 public:
  CompositePreconditioner(const std::shared_ptr<BaseLinearProgram> &lp) :
      BasePreconditioner(compose<Preconditioners...>(lp)) {}

  void initialize(const LinearProgram &prog) override { next()->initialize(prog); }

  void setObjective(const std::shared_ptr<Eigen::VectorXd> &c) override { next()->setObjective(c); }

  void setRHS(const std::shared_ptr<Eigen::VectorXd> &b_ub, const std::shared_ptr<Eigen::VectorXd> &b_eq) override {
    next()->setRHS(b_ub, b_eq);
  }

  void setBounds(const std::shared_ptr<Eigen::VectorXd> &lb, const std::shared_ptr<Eigen::VectorXd> &ub) override {
    next()->setBounds(lb, ub);
  }

  void addVar(double c, std::shared_ptr<Eigen::VectorXd> a_ub, std::shared_ptr<Eigen::VectorXd> a_eq, double lb,
              double ub) override { next()->addVar(c, a_ub, a_eq, lb, ub); }

  void addVars(std::shared_ptr<Eigen::VectorXd> c, std::shared_ptr<Eigen::SparseMatrix<double>> A_ub,
               std::shared_ptr<Eigen::SparseMatrix<double>> A_eq, std::shared_ptr<Eigen::VectorXd> lb,
               std::shared_ptr<Eigen::VectorXd> ub) override { next()->addVars(c, A_ub, A_eq, lb, ub); }

  void removeVar(size_t var) override { next()->removeVar(var); }

  void removeVars(size_t index, size_t n_remove) override { next()->removeVars(index, n_remove); }

  void removeBackVars(size_t n_remove) override { next()->removeBackVars(n_remove); }

  void addConstr_ub(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) override { next()->addConstr_ub(a, b); }

  void addConstrs_ub(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                     const std::shared_ptr<Eigen::VectorXd> &b) override { next()->addConstrs_ub(A, b); }

  void removeConstr_ub(size_t index) override { next()->removeConstr_ub(index); }

  void removeConstrs_ub(size_t index, size_t n_remove) override { next()->removeConstrs_ub(index, n_remove); }

  void addConstr_eq(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) override { next()->addConstr_eq(a, b); }

  void addConstrs_eq(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                     const std::shared_ptr<Eigen::VectorXd> &b) override { next()->addConstrs_eq(A, b); }

  void removeConstr_eq(size_t index) override { next()->removeConstr_eq(index); }

  void removeConstrs_eq(size_t index, size_t n_remove) override { next()->removeConstrs_eq(index, n_remove); }

  Eigen::VectorXd convertState(const Eigen::VectorXd &x) override { return next()->convertState(x); }

  Eigen::VectorXd revertState(const Eigen::VectorXd &x) override { return next()->revertState(x); }

 private:
  template<typename T>
  std::shared_ptr<BasePreconditioner> compose(const std::shared_ptr<BaseLinearProgram> &lp) {
    return std::make_shared<T>(lp);
  }

  template<typename T, typename T2, typename ...Ts>
  std::shared_ptr<BasePreconditioner> compose(const std::shared_ptr<BaseLinearProgram> &lp) {
    return std::make_shared<T>(compose<T2, Ts...>(lp));
  }
};
}}

#endif //STUKA_LP_COMPOSITE_PRECONDITIONER_H
