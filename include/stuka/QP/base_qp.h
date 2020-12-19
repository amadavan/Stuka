//
// Created by Avinash Madavan on 1/7/19.
//

#ifndef STUKA_QP_BASE_QP_H
#define STUKA_QP_BASE_QP_H

#include <memory>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../constants.h"
#include "qp.h"

namespace stuka { namespace QP {

class BaseQuadraticProgram {
 public:
  virtual ~BaseQuadraticProgram() = default;

  explicit BaseQuadraticProgram(const QuadraticProgram &prog) {};

  virtual void setObjective(const std::unique_ptr<Eigen::SparseMatrix<double>> &Q,
                            const std::unique_ptr<Eigen::VectorXd> &c) = 0;

  virtual void setRHS(const std::unique_ptr<Eigen::VectorXd> &b_ub, const std::unique_ptr<Eigen::VectorXd> &b_eq) = 0;

  virtual void setBounds(const std::unique_ptr<Eigen::VectorXd> &lb, const std::unique_ptr<Eigen::VectorXd> &ub) = 0;

  // TODO: take quadratic costs as well
  virtual void addVar(double c,
                      const std::unique_ptr<Eigen::VectorXd> &a_ub,
                      const std::unique_ptr<Eigen::VectorXd> &a_eq,
                      double lb,
                      double ub) = 0;

  virtual void addVars(const std::unique_ptr<Eigen::VectorXd> &c,
                       const std::unique_ptr<Eigen::SparseMatrix<double>> &A_ub,
                       const std::unique_ptr<Eigen::SparseMatrix<double>> &A_eq,
                       const std::unique_ptr<Eigen::VectorXd> &lb,
                       const std::unique_ptr<Eigen::VectorXd> &ub) = 0;

  virtual void removeVar(size_t var) = 0;

  virtual void removeVars(size_t index, size_t n_remove) = 0;

  virtual void removeBackVars(size_t n_remove) = 0;

  virtual void addConstr_ub(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) = 0;

  virtual void addConstrs_ub(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                             const std::unique_ptr<Eigen::VectorXd> &b) = 0;

  virtual void removeConstr_ub(size_t index) = 0;

  virtual void removeConstrs_ub(size_t index, size_t n_remove) = 0;

  virtual void addConstr_eq(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) = 0;

  virtual void addConstrs_eq(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                             const std::unique_ptr<Eigen::VectorXd> &b) = 0;

  virtual void removeConstr_eq(size_t index) = 0;

  virtual void removeConstrs_eq(size_t index, size_t n_remove) = 0;
};
}}

#endif //STUKA_QP_BASE_QP_H
