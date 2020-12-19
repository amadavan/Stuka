//
// Created by Avinash Madavan on 2019-03-02.
//

#ifndef STUKA_LP_SLACK_LP_H
#define STUKA_LP_SLACK_LP_H

#include <vector>

#include "base_lp.h"

namespace stuka { namespace LP {
class SlackLinearProgram : public BaseLinearProgram {
 public:
  void initialize(const LinearProgram &prog) override;

  void setObjective(const std::unique_ptr<Eigen::VectorXd> &c) override;

  void setRHS(const std::unique_ptr<Eigen::VectorXd> &b_ub, const std::unique_ptr<Eigen::VectorXd> &b_eq) override;

  void setBounds(const std::unique_ptr<Eigen::VectorXd> &lb, const std::unique_ptr<Eigen::VectorXd> &ub) override;

  void addVar(double c, const std::unique_ptr<Eigen::VectorXd> &a_ub, const std::unique_ptr<Eigen::VectorXd> &a_eq, double lb,
              double ub) override;

  void addVars(const std::unique_ptr<Eigen::VectorXd> &c, const std::unique_ptr<Eigen::SparseMatrix<double>> &A_ub,
               const std::unique_ptr<Eigen::SparseMatrix<double>> &A_eq, const std::unique_ptr<Eigen::VectorXd> &lb,
               const std::unique_ptr<Eigen::VectorXd> &ub) override;

  void removeVar(size_t var) override;

  void removeVars(size_t index, size_t n_remove) override;

  void removeBackVars(size_t n_remove) override;

  void addConstr_ub(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) override;

  void addConstrs_ub(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                     const std::unique_ptr<Eigen::VectorXd> &b) override;

  void removeConstr_ub(size_t index) override;

  void removeConstrs_ub(size_t index, size_t n_remove) override;

  void addConstr_eq(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) override;

  void addConstrs_eq(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                     const std::unique_ptr<Eigen::VectorXd> &b) override;

  void removeConstr_eq(size_t index) override;

  void removeConstrs_eq(size_t index, size_t n_remove) override;

  Eigen::VectorXd convertState(const Eigen::VectorXd &x) override;

  Eigen::VectorXd revertState(const Eigen::VectorXd &x) override;

 private:
  std::unique_ptr<Eigen::VectorXd> c_;
  std::unique_ptr<Eigen::SparseMatrix<double>> A_;
  std::unique_ptr<Eigen::VectorXd> b_;

  long n_dim_;
  long n_con_;
  long n_slack_;

  long n_split_;
  long n_const_;
  long n_add_;

  enum StateType {
    NORMAL,
    NEGATIVE,
    UNBOUNDED,
    CONSTANT,
    COMPACT
  };

  std::vector<StateType> x_type_;
  Eigen::VectorXd x_shift_;

  friend class MehrotraPC;
};
}}

#endif //STUKA_SLACK_LP_H
