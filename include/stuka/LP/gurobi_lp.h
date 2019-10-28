//
// Created by Avinash Madavan on 1/3/19.
//

#ifndef STUKA_LP_GUROBI_LP_H
#define STUKA_LP_GUROBI_LP_H

#include <algorithm>

#include <gurobi_c++.h>

#include "lp.h"
#include "base_lp.h"
#include "../constants.h"

namespace stuka { namespace LP {
class GurobiLinearProgram : public BaseLinearProgram {
 public:
  ~GurobiLinearProgram() override;

  GurobiLinearProgram();

  void setObjective(const std::shared_ptr<Eigen::VectorXd> &c) override;

  void setRHS(const std::shared_ptr<Eigen::VectorXd> &b_ub, const std::shared_ptr<Eigen::VectorXd> &b_eq) override;

  void setBounds(const std::shared_ptr<Eigen::VectorXd> &lb, const std::shared_ptr<Eigen::VectorXd> &ub) override;

  void addVar(double c,
              std::shared_ptr<Eigen::VectorXd> a_ub,
              std::shared_ptr<Eigen::VectorXd> a_eq,
              double lb,
              double ub) override;

  void addVars(std::shared_ptr<Eigen::VectorXd> c,
               std::shared_ptr<Eigen::SparseMatrix<double>> A_ub,
               std::shared_ptr<Eigen::SparseMatrix<double>> A_eq,
               std::shared_ptr<Eigen::VectorXd> lb,
               std::shared_ptr<Eigen::VectorXd> ub) override;

  void removeVar(size_t var) override;

  void removeVars(size_t index, size_t n_remove) override;

  void removeBackVars(size_t n_remove) override;

  void addConstr_ub(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) override;

  void addConstrs_ub(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                     const std::shared_ptr<Eigen::VectorXd> &b) override;

  void removeConstr_ub(size_t index) override;

  void removeConstrs_ub(size_t index, size_t n_remove) override;

  void addConstr_eq(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) override;

  void addConstrs_eq(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                     const std::shared_ptr<Eigen::VectorXd> &b) override;

  void removeConstr_eq(size_t index) override;

  void removeConstrs_eq(size_t index, size_t n_remove) override;

  Eigen::VectorXd convertState(const Eigen::VectorXd &x) override;

  Eigen::VectorXd revertState(const Eigen::VectorXd &x) override;

  void initialize(const LinearProgram &prog) override;

 private:
  GRBEnv env_;
  GRBModel model_;

  GRBVar *vars_;
  GRBConstr *eqconstr_;
  GRBConstr *ubconstr_;

  long n_dim_;
  long n_con_ub_;
  long n_con_eq_;

  long n_alloc_;                // Number of allocated variables (to prevent unnecessary copies)
  long n_alloc_ub_;             // Number of allocated ub constraints (to prevent unnecessary copies)
  long n_alloc_eq_;             // Number of allocated eq constraints (to prevent unnecessary copies)

  friend class GurobiSolver;
};
}}

#endif //STUKA_LP_GUROBI_LP_H
