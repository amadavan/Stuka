//
// Created by Avinash on 8/19/2020.
//

#ifndef STUKA_INCLUDE_STUKA_LP_LAZY_LP_H_
#define STUKA_INCLUDE_STUKA_LP_LAZY_LP_H_

#include <iostream>

#include "base_lp.h"

namespace stuka { namespace LP {

class LazyLinearProgram : public BaseLinearProgram {
 public:
  LazyLinearProgram();

  ~LazyLinearProgram() override = default;

  void initialize(const LinearProgram &prog) override;

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

  void setLP(BaseLinearProgram *prog);

  bool isActive(size_t index);

  void setConstrActive(size_t index);

  void setConstrsActive(size_t index, size_t n_active);

  void setConstrsActive(Eigen::Matrix<bool, Eigen::Dynamic, 1> activity);

  void setConstrsActive(std::vector<size_t> activity);

  bool isViolated(const Eigen::VectorXd &x);

  bool addViolatedConstraints(const Eigen::VectorXd &x);

  bool addRandomConstraint();

  Eigen::VectorXd getDualUB(const Eigen::VectorXd &active_dual);

 private:
  std::shared_ptr<Eigen::SparseMatrix<double>> A_ub_;
  std::shared_ptr<Eigen::VectorXd> b_ub_;

  size_t n_ub_;
  size_t n_active_;
  Eigen::Matrix<bool, Eigen::Dynamic, 1> ub_activity_;
  Eigen::Matrix<size_t, Eigen::Dynamic, 1> ub_index_;

  BaseLinearProgram *prog_;

};

}}

#endif //STUKA_INCLUDE_STUKA_LP_LAZY_LP_H_
