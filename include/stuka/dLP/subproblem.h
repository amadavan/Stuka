//
// Created by Avinash Madavan on 1/15/19.
//

#ifndef STUKA_dLP_SUBPROBLEM_H
#define STUKA_dLP_SUBPROBLEM_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../LP/lp.h"

namespace stuka { namespace dLP {
struct Subproblem {
  std::unique_ptr<Eigen::VectorXd> c;
  std::unique_ptr<Eigen::SparseMatrix<double>> A_ub;
  std::unique_ptr<Eigen::VectorXd> b_ub;
  std::unique_ptr<Eigen::SparseMatrix<double>> C_ub;
  std::unique_ptr<Eigen::SparseMatrix<double>> A_eq;
  std::unique_ptr<Eigen::VectorXd> b_eq;
  std::unique_ptr<Eigen::SparseMatrix<double>> C_eq;
  std::unique_ptr<Eigen::VectorXd> lb;
  std::unique_ptr<Eigen::VectorXd> ub;
};
}}

#endif //STUKA_dLP_SUBPROBLEM_H
