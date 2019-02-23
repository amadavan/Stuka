//
// Created by Avinash Madavan on 7/16/18.
//

#include "simple_lp.h"

stuka::LP::LinearProgram stuka::example::SimpleLP::gen() {
  stuka::LP::LinearProgram prog;

  prog.c = std::make_shared<Eigen::VectorXd>(4);
  Eigen::MatrixXd A_ub = Eigen::MatrixXd(2, 4);
  prog.b_ub = std::make_shared<Eigen::VectorXd>(2);
  Eigen::MatrixXd A_eq = Eigen::MatrixXd(1, 4);
  prog.b_eq = std::make_shared<Eigen::VectorXd>(1);
  prog.lb = std::make_shared<Eigen::VectorXd>(4);
  prog.ub = std::make_shared<Eigen::VectorXd>(4);

  *prog.c << 2, 1, 0, 0;
  A_ub << -1, -1, 0, 0, 1, -2, 0, 0;
  *prog.b_ub << -2, 4;
  A_eq << -1, 1, 1, 0;
  *prog.b_eq << 1;
  *prog.lb << -stuka::INF, 0, 0, 0;
  *prog.ub << stuka::INF, stuka::INF, stuka::INF, 0;

  prog.A_ub = std::make_shared<Eigen::SparseMatrix<double>>(A_ub.sparseView());
  prog.A_eq = std::make_shared<Eigen::SparseMatrix<double>>(A_eq.sparseView());

  return prog;
}
