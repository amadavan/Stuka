//
// Created by Avinash Madavan on 10/8/19.
//

#include <stuka/LP/lp.h>

void stuka::LP::LinearProgram::ConstraintReductionBounds() {
  if (!b_ub || !lb || !ub) return;

  Eigen::SparseMatrix<double> zero(A_ub->rows(), A_ub->cols());
  zero.reserve(A_ub->nonZeros());
  for (size_t i = 0; i < A_ub->cols(); ++i) {
    zero.startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub, i); it; ++it)
      zero.insertBack(it.row(), i) = 0;
  }

  Eigen::VectorXd U = A_ub->cwiseMax(zero) * *ub + A_ub->cwiseMin(zero) * *lb;
  Eigen::Matrix<bool, Eigen::Dynamic, 1> redundancy = (*b_ub - U).array() > 0;

  size_t n_dim = c->size();
  size_t n_ub = b_ub->size();
  size_t n_redundant = redundancy.count();

  std::shared_ptr<Eigen::SparseMatrix<double>>
      A_ub_reduced = std::make_shared<Eigen::SparseMatrix<double>>(n_ub - n_redundant, n_dim);
  std::shared_ptr<Eigen::VectorXd> b_ub_reduced = std::make_shared<Eigen::VectorXd>(n_ub - n_redundant);
  A_ub_reduced->reserve(A_ub->nonZeros());

  // Count all previous active constraints at each index. This will be used in the next step to construct the complete
  // sets of active and inactive constraints
  Eigen::VectorXi redundant_count(n_ub);
  Eigen::VectorXi active_count(n_ub);

  int count_redundant = 0, count_active = 0;
  for (size_t i = 0; i < n_ub; ++i) {
    redundant_count.coeffRef(i) = count_redundant;
    active_count.coeffRef(i) = count_active;
    if (redundancy.coeff(i)) ++count_redundant;
    else ++count_active;
  }

  for (size_t i = 0; i < n_dim; ++i) {
    A_ub_reduced->startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub, i); it; ++it)
      if (!redundancy.coeff(it.row()))
        A_ub_reduced->insertBack(active_count.coeff(it.row()), i) = it.value();
  }

  for (size_t i = 0; i < n_ub; ++i)
    if (!redundancy.coeff(i))
      b_ub_reduced->coeffRef(active_count.coeff(i)) = b_ub->coeff(i);

  A_ub = A_ub_reduced;
  b_ub = b_ub_reduced;
}
