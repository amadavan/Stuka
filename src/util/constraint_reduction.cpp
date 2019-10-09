//
// Created by Avinash Madavan on 10/8/19.
//

#include <stuka/LP/constraint_reduction.h>

void stuka::LP::ConstraintReduction::bounds(const LinearProgram &prog) {
  if (!prog.b_ub) return;

  Eigen::VectorXd U = prog.A_ub->cwiseMax(0) * prog.ub + prog.A_ub->cwiseMin(0) * prog.lb;
  Eigen::Matrix<bool, Eigen::Dynamic, 1> activity = (*prog.b_ub - U).array() > 0;

  size_t n_dim = prog.c->size();
  size_t n_ub = prog.b_ub->size();
  size_t n_redundant = activity.count();

  std::shared_ptr<Eigen::SparseMatrix<double>> A_ub = std::make_shared<Eigen::SparseMatrix<double>>(n_ub - n_redundant, n_dim);
  std::shared_ptr<Eigen::VectorXd> b_ub = std::make_shared<Eigen::VectorXd>(n_ub - n_redundant);

  // Count all previous active constraints at each index. This will be used in the next step to construct the complete
  // sets of active and inactive constraints
  Eigen::VectorXi active_count(n_ub);
  Eigen::VectorXi inactive_count(n_ub);

  int count_active = 0;
  for (size_t i = 0; i < n_ub; ++i) {
    active_count.coeffRef(i) = count_active;
    if (activity.coeff(i)) ++count_active;
  }

  for (size_t i = 0; i < n_dim; ++i) {
    A_ub->startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_ub, i); it; ++it)
      if (activity.coeff(it.row()))
        A_ub.insertBack(active_count.coeff(it.row()), i) = it.value();
  }

  for (size_t i = 0; i < n_ub; ++i)
    if (activity.coeff(i))
      b_ub.coeffRef(active_count.coeff(i)) = prog.b_ub->coeff(i);

  prog.A_ub = A_ub;
  prog.b_ub = b_ub;
}
