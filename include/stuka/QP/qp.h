//
// Created by Avinash Madavan on 1/7/19.
//

#ifndef STUKA_QP_QP_H
#define STUKA_QP_QP_H

#include <memory>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace stuka { namespace QP {
struct QuadraticProgram {
  std::shared_ptr<Eigen::SparseMatrix<double>> Q;
  std::shared_ptr<Eigen::VectorXd> c;
  std::shared_ptr<Eigen::SparseMatrix<double>> A_ub;
  std::shared_ptr<Eigen::VectorXd> b_ub;
  std::shared_ptr<Eigen::SparseMatrix<double>> A_eq;
  std::shared_ptr<Eigen::VectorXd> b_eq;
  std::shared_ptr<Eigen::VectorXd> lb;
  std::shared_ptr<Eigen::VectorXd> ub;

  // Hack for cleaner python implementation
  virtual std::string name() { return "QuadraticProgram"; }
};
}}

#endif //STUKA_QP_QP_H
