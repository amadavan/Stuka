//
// Created by Avinash Madavan on 1/19/20.
//

#include <stuka/util/constraint_reduction.h>

std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd>
stuka::util::LinearInequalityConstraintReduction::bounds(Eigen::VectorXd lb, Eigen::VectorXd ub) {

  Eigen::SparseMatrix<double> zero(A_.rows(), A_.cols());
  zero.reserve(A_.nonZeros());
  for (size_t i = 0; i < A_.cols(); ++i) {
    zero.startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(A_, i); it; ++it)
      zero.insertBack(it.row(), i) = 0;
  }

  Eigen::VectorXd U = A_.cwiseMax(zero) * ub + A_.cwiseMin(zero) * lb;
  Eigen::Matrix<bool, Eigen::Dynamic, 1> redundancy = (b_ - U).array() > 0;

  size_t n_dim = A_.rows();
  size_t n_ub = b_.size();
  size_t n_redundant = redundancy.count();

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

  Eigen::SparseMatrix<double> A_reduced = Eigen::SparseMatrix<double>(n_ub - n_redundant, n_dim);
  Eigen::VectorXd b_reduced = Eigen::VectorXd::Zero(n_ub - n_redundant);

  A_reduced.reserve(A_.nonZeros());
  for (size_t i = 0; i < n_dim; ++i) {
    A_reduced.startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(A_, i); it; ++it)
      if (!redundancy.coeff(it.row()))
        A_reduced.insertBack(active_count.coeff(it.row()), i) = it.value();
  }
  A_reduced.finalize();

  for (size_t i = 0; i < n_ub; ++i)
    if (!redundancy.coeff(i))
      b_reduced.coeffRef(active_count.coeff(i)) = b_.coeff(i);

  return std::make_pair(A_reduced, b_reduced);
}

std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd>
stuka::util::LinearInequalityConstraintReduction::stojkovic_stanimirovic(Eigen::VectorXd c) {
  Eigen::SparseMatrix<double> D = Eigen::SparseMatrix<double>(A_);

  for (size_t i = 0; i < c.size(); ++i)
    for (Eigen::SparseMatrix<double>::InnerIterator it(D, i); it; ++it)
      D.coeffRef(i, it.col()) /= c.coeff(i) * b_.coeff(it.col());

  // Need a proper way to compare two rows of a sparse matrix
  for (size_t i = 0; i < D.rows() - 1; ++i) {
    for (size_t j = i+1; j < D.rows(); ++j) {
      for (size_t k = 0; k < D.cols(); ++k) {
        // TODO: complete this.
      }
    }
  }
}
