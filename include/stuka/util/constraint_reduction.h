//
// Created by Avinash Madavan on 1/2/20.
//

#ifndef STUKA_UTIL_CONSTRAINT_REDUCTION_H_
#define STUKA_UTIL_CONSTRAINT_REDUCTION_H_

#include <Eigen/SparseCore>

#include "../LP/lp.h"

namespace stuka { namespace util {

class LinearInequalityConstraintReduction {
 public:
  LinearInequalityConstraintReduction(const Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b) : A_(A), b_(b) {}

  ~LinearInequalityConstraintReduction() {}

  std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> bounds(Eigen::VectorXd lb, Eigen::VectorXd ub);

  std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> stojkovic_stanimirovic(Eigen::VectorXd c);

  std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> paulraj();

 private:
  const Eigen::SparseMatrix<double> &A_;
  const Eigen::VectorXd &b_;
};

}}

#endif //STUKA_UTIL_CONSTRAINT_REDUCTION_H_
