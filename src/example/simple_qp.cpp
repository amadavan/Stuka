//
// Created by Avinash Madavan on 7/18/18.
//

#include "simple_qp.h"

stuka::QP::QuadraticProgram stuka::example::SimpleQP::gen() {
  stuka::QP::QuadraticProgram prog;

  prog.Q = std::make_shared<Eigen::SparseMatrix<double>>(2, 2);
  prog.Q->coeffRef(0, 0) = 1;
  prog.Q->coeffRef(0, 1) = -1;
  prog.Q->coeffRef(1, 0) = -1;
  prog.Q->coeffRef(1, 1) = 2;

  prog.c = std::make_shared<Eigen::VectorXd>(2);
  *prog.c << -2, -6;

  prog.A_ub = std::make_shared<Eigen::SparseMatrix<double>>(3, 2);
  prog.A_ub->coeffRef(0, 0) = 1;
  prog.A_ub->coeffRef(1, 0) = -1;
  prog.A_ub->coeffRef(2, 0) = 2;
  prog.A_ub->coeffRef(0, 1) = 1;
  prog.A_ub->coeffRef(1, 1) = 2;
  prog.A_ub->coeffRef(2, 1) = 1;

  prog.b_ub = std::make_shared<Eigen::VectorXd>(3);
  *prog.b_ub << 2, 2, 3;

  prog.lb = std::make_shared<Eigen::VectorXd>(2);
  *prog.lb << 0, 0;

  prog.ub = std::make_shared<Eigen::VectorXd>(2);
  *prog.ub << stuka::INF, stuka::INF;

  return prog;
}

std::string stuka::example::SimpleQP::name() {
  return std::string("MATLAB Example Quadratic Program");
}
