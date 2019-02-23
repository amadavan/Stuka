//
// Created by Avinash Madavan on 7/16/18.
//

#include "borrelli_729.h"
#include <iostream>

stuka::dLP::DecomposedLinearProgram stuka::example::Borrelli729::gen() {
  typedef Eigen::Triplet<double> T;

  dLP::DecomposedLinearProgram dlp(1);

  // Master problem
  dlp.c[1] = std::make_shared<Eigen::VectorXd>(2);

  *dlp.c[1] << 0, 0;

  dlp.A_ub[1] = nullptr;
  dlp.b_ub[1] = nullptr;
  dlp.A_eq[1] = nullptr;
  dlp.b_eq[1] = nullptr;

  dlp.lb[1] = std::make_shared<Eigen::VectorXd>(2);
  *dlp.lb[1] << -2.5, -2.5;

  dlp.ub[1] = std::make_shared<Eigen::VectorXd>(2);
  *dlp.ub[1] << 2.5, 2.5;
  // Subproblem
  dlp.c[0] = std::make_shared<Eigen::VectorXd>(4);
  *dlp.c[0] << 1, 1, 0, 0;

  dlp.A_ub[0] = std::make_shared<Eigen::SparseMatrix<double>>(8, 4);
  std::vector<T> A_ub1_triplet(16);
  A_ub1_triplet.emplace_back(T(0, 0, -1.));
  A_ub1_triplet.emplace_back(T(1, 0, -1.));
  A_ub1_triplet.emplace_back(T(2, 0, -1.));
  A_ub1_triplet.emplace_back(T(3, 0, -1.));
  A_ub1_triplet.emplace_back(T(4, 1, -1.));
  A_ub1_triplet.emplace_back(T(5, 1, -1.));
  A_ub1_triplet.emplace_back(T(6, 1, -1.));
  A_ub1_triplet.emplace_back(T(7, 1, -1.));
  A_ub1_triplet.emplace_back(T(1, 2, -1.));
  A_ub1_triplet.emplace_back(T(3, 2, 1.));
  A_ub1_triplet.emplace_back(T(4, 2, -1.));
  A_ub1_triplet.emplace_back(T(5, 2, -1.));
  A_ub1_triplet.emplace_back(T(6, 2, 1.));
  A_ub1_triplet.emplace_back(T(7, 2, 1.));
  A_ub1_triplet.emplace_back(T(5, 3, -1.));
  A_ub1_triplet.emplace_back(T(7, 3, 1.));
  dlp.A_ub[0]->setFromTriplets(A_ub1_triplet.begin(), A_ub1_triplet.end());

  dlp.b_ub[0] = std::make_shared<Eigen::VectorXd>(8);
  dlp.b_ub[0]->setZero();

  dlp.C_ub[0] = std::make_shared<Eigen::SparseMatrix<double>>(8, 2);
  std::vector<T> C_ub1_triplet(12);
  C_ub1_triplet.emplace_back(T(0, 0, -1.));
  C_ub1_triplet.emplace_back(T(2, 0, 1.));
  C_ub1_triplet.emplace_back(T(4, 0, -1.));
  C_ub1_triplet.emplace_back(T(6, 0, 1.));
  C_ub1_triplet.emplace_back(T(0, 1, -1.));
  C_ub1_triplet.emplace_back(T(1, 1, -1.));
  C_ub1_triplet.emplace_back(T(2, 1, 1.));
  C_ub1_triplet.emplace_back(T(3, 1, 1.));
  C_ub1_triplet.emplace_back(T(4, 1, -2.));
  C_ub1_triplet.emplace_back(T(5, 1, -1.));
  C_ub1_triplet.emplace_back(T(6, 1, 2.));
  C_ub1_triplet.emplace_back(T(7, 1, 1.));
  dlp.C_ub[0]->setFromTriplets(C_ub1_triplet.begin(), C_ub1_triplet.end());

  dlp.A_eq[0] = nullptr;
  dlp.b_eq[0] = nullptr;
  dlp.C_eq[0] = nullptr;

  dlp.lb[0] = std::make_shared<Eigen::VectorXd>(4);
  *dlp.lb[0] << -stuka::INF, -stuka::INF, -1, -1;

  dlp.ub[0] = std::make_shared<Eigen::VectorXd>(4);
  *dlp.ub[0] << stuka::INF, stuka::INF, 1, 1;

  return dlp;
}

stuka::LP::LinearProgram stuka::example::Borrelli729::full() {
  typedef Eigen::Triplet<double> T;

  LP::LinearProgram lp;

  lp.c = std::make_shared<Eigen::VectorXd>(6);
  *lp.c << 1, 1, 0, 0, 0, 0;

  Eigen::SparseMatrix<double> A_ub1 = Eigen::SparseMatrix<double>(8, 4);
  Eigen::SparseMatrix<double> C_ub1 = Eigen::SparseMatrix<double>(8, 2);

  std::vector<T> A_ub1_triplet(16);
  A_ub1_triplet.emplace_back(T(0, 0, -1.));
  A_ub1_triplet.emplace_back(T(1, 0, -1.));
  A_ub1_triplet.emplace_back(T(2, 0, -1.));
  A_ub1_triplet.emplace_back(T(3, 0, -1.));
  A_ub1_triplet.emplace_back(T(4, 1, -1.));
  A_ub1_triplet.emplace_back(T(5, 1, -1.));
  A_ub1_triplet.emplace_back(T(6, 1, -1.));
  A_ub1_triplet.emplace_back(T(7, 1, -1.));
  A_ub1_triplet.emplace_back(T(1, 2, -1.));
  A_ub1_triplet.emplace_back(T(3, 2, 1.));
  A_ub1_triplet.emplace_back(T(4, 2, -1.));
  A_ub1_triplet.emplace_back(T(5, 2, -1.));
  A_ub1_triplet.emplace_back(T(6, 2, 1.));
  A_ub1_triplet.emplace_back(T(7, 2, 1.));
  A_ub1_triplet.emplace_back(T(5, 3, -1.));
  A_ub1_triplet.emplace_back(T(7, 3, 1.));
  A_ub1.setFromTriplets(A_ub1_triplet.begin(), A_ub1_triplet.end());

  std::vector<T> C_ub1_triplet(12);
  C_ub1_triplet.emplace_back(T(0, 0, -1.));
  C_ub1_triplet.emplace_back(T(2, 0, 1.));
  C_ub1_triplet.emplace_back(T(4, 0, -1.));
  C_ub1_triplet.emplace_back(T(6, 0, 1.));
  C_ub1_triplet.emplace_back(T(0, 1, -1.));
  C_ub1_triplet.emplace_back(T(1, 1, -1.));
  C_ub1_triplet.emplace_back(T(2, 1, 1.));
  C_ub1_triplet.emplace_back(T(3, 1, 1.));
  C_ub1_triplet.emplace_back(T(4, 1, -2.));
  C_ub1_triplet.emplace_back(T(5, 1, -1.));
  C_ub1_triplet.emplace_back(T(6, 1, 2.));
  C_ub1_triplet.emplace_back(T(7, 1, 1.));
  C_ub1.setFromTriplets(C_ub1_triplet.begin(), C_ub1_triplet.end());

  lp.A_ub = std::make_shared<Eigen::SparseMatrix<double>>(8, 6);

  lp.A_ub->setZero();
  for (int i = 0; i < 4; ++i) {
    lp.A_ub->startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(A_ub1, i); it; ++it)
      lp.A_ub->insertBack(it.row(), i) = it.value();
  }
  for (int i = 0; i < 2; ++i) {
    lp.A_ub->startVec(i + 4);
    for (Eigen::SparseMatrix<double>::InnerIterator it(C_ub1, i); it; ++it)
      lp.A_ub->insertBack(it.row(), i + 4) = it.value();
  }
  lp.A_ub->finalize();

  lp.b_ub = std::make_shared<Eigen::VectorXd>(8);
  lp.b_ub->setZero();

  lp.A_eq = nullptr;
  lp.b_eq = nullptr;

  lp.lb = std::make_shared<Eigen::VectorXd>(6);
  *lp.lb << -stuka::INF, -stuka::INF, -1, -1, -2.5, -2.5;

  lp.ub = std::make_shared<Eigen::VectorXd>(6);
  *lp.ub << stuka::INF, stuka::INF, 1, 1, 2.5, 2.5;

  return lp;
}
