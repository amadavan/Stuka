//
// Created by Avinash Madavan on 2019-07-06.
//

#ifndef STUKA_LP_BOUNDS_ZERO_H
#define STUKA_LP_BOUNDS_ZERO_H

#include "base_preconditioner.h"

namespace stuka { namespace LP {
class BoundsZero : public BasePreconditioner {
 public:
  BoundsZero(const std::shared_ptr<BaseLinearProgram> &next);

  void initialize(const LinearProgram &prog) override;

  void setObjective(const std::shared_ptr<Eigen::VectorXd> &c) override;

  void setRHS(const std::shared_ptr<Eigen::VectorXd> &b_ub, const std::shared_ptr<Eigen::VectorXd> &b_eq) override;

  void setBounds(const std::shared_ptr<Eigen::VectorXd> &lb, const std::shared_ptr<Eigen::VectorXd> &ub) override;

  void addVar(double c, std::shared_ptr<Eigen::VectorXd> a_ub, std::shared_ptr<Eigen::VectorXd> a_eq, double lb,
              double ub) override;

  void addVars(std::shared_ptr<Eigen::VectorXd> c, std::shared_ptr<Eigen::SparseMatrix<double>> A_ub,
               std::shared_ptr<Eigen::SparseMatrix<double>> A_eq, std::shared_ptr<Eigen::VectorXd> lb,
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

 private:
  size_t n_dim_;
  size_t n_ub_;
  size_t n_split_;
  size_t n_bounds_;
  size_t n_slack_;

  enum StateType {
    NORMAL,
    NEGATIVE,
    UNBOUNDED,
    COMPACT
  };

  Eigen::VectorXd x_shift_;
  std::vector<StateType> x_type_;
};
}}

#endif //STUKA_LP_BOUNDS_ZERO_H
