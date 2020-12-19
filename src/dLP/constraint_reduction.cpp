//
// Created by Avinash Madavan on 10/8/19.
//

#include <iostream>

#include <stuka/dLP/decomposed_lp.h>
#include <stuka/LP/lp.h>

void stuka::dLP::DecomposedLinearProgram::ConstraintReductionBounds() {
  if (!lb.back() || !ub.back()) return;

  size_t n_dim0 = c.back()->size();

  LP::LinearProgram master;
  master.c = util::DenseOps::unique_copy(c.back());
  master.A_ub = util::SparseOps::unique_copy(A_ub.back());
  master.b_ub = util::DenseOps::unique_copy(b_ub.back());
  master.A_eq = util::SparseOps::unique_copy(A_eq.back());
  master.b_eq = util::DenseOps::unique_copy(b_eq.back());
  master.lb = util::DenseOps::unique_copy(lb.back());
  master.ub = util::DenseOps::unique_copy(ub.back());

  master.reduceConstraints(ConstraintReductionMethods::BOUNDS);

  A_ub.back() = std::move(master.A_ub);
  b_ub.back() = std::move(master.b_ub);

  size_t n_sub = c.size() - 1;
  for (size_t i = 0; i < n_sub; ++i) {
    if (!lb[i] || !ub[i]) continue;

    Eigen::SparseMatrix<double> zero_A(A_ub[i]->rows(), A_ub[i]->cols());
    zero_A.reserve(A_ub[i]->nonZeros());
    for (size_t j = 0; j < A_ub[i]->cols(); ++j) {
      zero_A.startVec(j);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub[i], j); it; ++it)
        zero_A.insertBack(it.row(), j) = 0;
    }

    Eigen::SparseMatrix<double> zero_C(C_ub[i]->rows(), C_ub[i]->cols());
    zero_C.reserve(C_ub[i]->nonZeros());
    for (size_t j = 0; j < C_ub[j + 1]->cols(); ++j) {
      zero_C.startVec(j);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*C_ub[j + 1], j); it; ++it)
        zero_C.insertBack(it.row(), j) = 0;
    }

    Eigen::VectorXd Ui = A_ub[i]->cwiseMax(zero_A) * *ub[i] + A_ub[i]->cwiseMin(zero_A) * *lb[i];
    Eigen::VectorXd U0 = C_ub[i]->cwiseMax(zero_C) * *ub.back() + C_ub[i]->cwiseMin(zero_C) * *lb.back();
    Eigen::Array<bool, Eigen::Dynamic, 1> redundancy = (*(b_ub[i]) - U0 - Ui).array() > 0;

    size_t n_dim = c[i]->size();
    size_t n_ub = b_ub[i]->size();
    size_t n_redundant = redundancy.count();

    std::unique_ptr<Eigen::SparseMatrix<double>>
        A_ub_reduced = std::make_unique<Eigen::SparseMatrix<double>>(n_ub - n_redundant, n_dim);
    std::unique_ptr<Eigen::VectorXd> b_ub_reduced = std::make_unique<Eigen::VectorXd>(n_ub - n_redundant);
    std::unique_ptr<Eigen::SparseMatrix<double>>
        C_ub_reduced = std::make_unique<Eigen::SparseMatrix<double>>(n_ub - n_redundant, n_dim0);

    A_ub_reduced->reserve(A_ub[i]->nonZeros());
    C_ub_reduced->reserve(C_ub[i]->nonZeros());

    // Count all previous active constraints at each index. This will be used in the next step to construct the complete
    // sets of active and inactive constraints
    Eigen::VectorXi redundant_count(n_ub);
    Eigen::VectorXi active_count(n_ub);

    int count_redundant = 0, count_active = 0;
    for (size_t j = 0; j < n_ub; ++j) {
      redundant_count.coeffRef(j) = count_redundant;
      active_count.coeffRef(j) = count_active;
      if (redundancy.coeff(j)) ++count_redundant;
      else ++count_active;
    }

    for (size_t j = 0; j < n_dim; ++j) {
      A_ub_reduced->startVec(j);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub[i], j); it; ++it)
        if (!redundancy.coeff(it.row()))
          A_ub_reduced->insertBack(active_count.coeff(it.row()), j) = it.value();
    }

    for (size_t j = 0; j < n_ub; ++j)
      if (!redundancy.coeff(j))
        b_ub_reduced->coeffRef(active_count.coeff(j)) = b_ub[i]->coeff(j);

    for (size_t j = 0; j < n_dim0; ++j) {
      C_ub_reduced->startVec(j);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*C_ub[i], j); it; ++it)
        if (!redundancy.coeff(it.row()))
          C_ub_reduced->insertBack(active_count.coeff(it.row()), j) = it.value();
    }

    A_ub[i] = std::move(A_ub_reduced);
    b_ub[i] = std::move(b_ub_reduced);
    C_ub[i] = std::move(C_ub_reduced);

    std::cout << n_redundant << std::endl;
  }

  return;
}

