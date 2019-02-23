//
// Created by Avinash Madavan on 1/15/19.
//

#include <stuka/dLP/cre_subproblem.h>

stuka::dLP::CRESubproblem::CRESubproblem(const stuka::dLP::Subproblem sub) : sub_(sub) {

  n_dim_ = sub_.c->size();
  n_dim_master_ = (sub_.C_ub) ? sub_.C_ub->cols() : (sub_.C_eq) ? sub_.C_eq->cols() : 0;

  n_bounds_ = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
    if (sub_.lb && sub_.lb->coeff(i) > -INF) n_bounds_++;
    if (sub_.ub && sub_.ub->coeff(i) < INF) n_bounds_++;
  }

  LP::LinearProgram lp;
  lp.c = sub_.c;
  lp.A_ub = sub_.A_ub;
  lp.b_ub = sub_.b_ub;
  lp.A_eq = sub_.A_eq;
  lp.b_eq = sub_.b_eq;
  lp.lb = sub_.lb;
  lp.ub = sub_.ub;

  // TODO: use solver selected by options
  solver_ = util::createSolver(lp);
}

stuka::dLP::CriticalRegion stuka::dLP::CRESubproblem::computeCriticalRegion(const Eigen::VectorXd &x) const {

  solver_->getLP().setRHS((sub_.b_ub) ? std::make_shared<Eigen::VectorXd>(*sub_.b_ub - *sub_.C_ub * x) : nullptr,
                          (sub_.b_eq) ? std::make_shared<Eigen::VectorXd>(*sub_.b_eq - *sub_.C_eq * x) : nullptr);

  // Solve the subproblem at the given point
  OptimizeState res;
  try {
    res = solver_->solve();
  } catch (std::exception &e) {
//    std::cout << "Failed to solve subproblem: " << e.getMessage() << std::endl;
//    std::cout << e.getErrorCode() << std::endl;
    throw std::runtime_error("computeCriticalRegion: unable to solve subproblem");
  }
  if (res.status != 2) {
    std::cout << "Failed to solve subproblem: " << res.status << std::endl;
    throw std::runtime_error("computeCriticalRegion: unable to solve subproblem");
  }

  Eigen::Matrix<bool, Eigen::Dynamic, 1> ub_activity = res.dual_ub.array() > 0;
  Eigen::Matrix<bool, Eigen::Dynamic, 1> xlb_activity = res.dual_x_lb.array() > 0;
  Eigen::Matrix<bool, Eigen::Dynamic, 1> xub_activity = res.dual_x_ub.array() > 0;

  long n_con_ub = (sub_.b_ub) ? sub_.b_ub->size() : 0;

  // Count all previous active constraints at each index. This will be used in the next step to construct the complete
  // sets of active and inactive constraints
  Eigen::VectorXi active_count_ub(n_con_ub);
  Eigen::VectorXi inactive_count_ub(n_con_ub);

  int count_active = 0, count_inactive = 0;
  for (size_t i = 0; i < n_con_ub; ++i) {
    active_count_ub.coeffRef(i) = count_active;
    inactive_count_ub.coeffRef(i) = count_inactive;
    if (ub_activity.coeff(i)) ++count_active;
    else ++count_inactive;
  }

  // Count of active constraints
  long n_eq = (sub_.b_eq) ? sub_.b_eq->size() : 0;
  long n_ub = (sub_.b_ub) ? sub_.b_ub->size() : 0;
  long n_ub_active = ub_activity.count();
  long n_ub_inactive = n_ub - ub_activity.count();
  long n_xlb = xlb_activity.count();
  long n_xub = xub_activity.count();

  long n_active = n_eq + n_ub_active + n_xlb + n_xub;
  long n_inactive = n_ub_inactive + n_bounds_ - n_xlb - n_xub;

  // Construct set of active constraints -------------------------------------------------------------------------------
  size_t nnz;

  // Active A matrix
  Eigen::SparseMatrix<double> A_active(n_active, n_dim_);
  nnz = n_xlb + n_xub;
  if (sub_.A_eq) nnz += sub_.A_eq->nonZeros();
  if (sub_.A_ub) nnz += sub_.A_ub->nonZeros();
  A_active.reserve(nnz);

  size_t n_xbound = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
    A_active.startVec(i);
    if (n_eq > 0)
      for (Eigen::SparseMatrix<double>::InnerIterator it(*sub_.A_eq, i); it; ++it)
        A_active.insertBack(it.row(), i) = it.value();
    if (n_ub_active > 0)
      for (Eigen::SparseMatrix<double>::InnerIterator it(*sub_.A_ub, i); it; ++it)
        if (ub_activity.coeff(it.row()))
          A_active.insertBack(n_eq + active_count_ub.coeff(it.row()), i) = it.value();
    if (xlb_activity[i]) {
      A_active.insertBack(n_eq + n_ub_active + n_xbound, i) = -1;
      n_xbound++;
    } else if (xub_activity[i]) {
      A_active.insertBack(n_eq + n_ub_active + n_xbound, i) = 1;
      n_xbound++;
    }
  }
  A_active.finalize();

  // Active b vector
  Eigen::VectorXd b_active(n_active);
  if (sub_.b_eq) b_active.head(n_eq) = *sub_.b_eq;
  for (size_t i = 0; i < n_ub; ++i)
    if (ub_activity.coeff(i))
      b_active.coeffRef(n_eq + active_count_ub.coeff(i)) = sub_.b_ub->coeff(i);
  n_xbound = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
    if (sub_.lb && xlb_activity[i]) {
      b_active.coeffRef(n_eq + n_ub_active + n_xbound) = -sub_.lb->coeff(i);
      n_xbound++;
    } else if (sub_.ub && xub_activity[i]) {
      b_active.coeffRef(n_eq + n_ub_active + n_xbound) = sub_.ub->coeff(i);
      n_xbound++;
    }
  }

  // Active C matrix
  Eigen::SparseMatrix<double> C_active(n_active, n_dim_master_);
  nnz = 0;
  if (sub_.C_eq) nnz += sub_.C_eq->nonZeros();
  if (sub_.C_ub) nnz += sub_.C_ub->nonZeros();
  C_active.reserve(nnz);
  for (size_t i = 0; i < n_dim_master_; ++i) {
    C_active.startVec(i);
    if (n_eq > 0)
      for (Eigen::SparseMatrix<double>::InnerIterator it(*sub_.C_eq, i); it; ++it)
        C_active.insertBack(it.row(), i) = it.value();
    if (n_ub_active > 0)
      for (Eigen::SparseMatrix<double>::InnerIterator it(*sub_.C_ub, i); it; ++it)
        if (ub_activity.coeff(it.row()))
          C_active.insertBack(n_eq + active_count_ub.coeff(it.row()), i) = it.value();
  }
  C_active.finalize();

  // Construct set of inactive constraints -----------------------------------------------------------------------------
  // Inactive A matrix
  Eigen::SparseMatrix<double, 0, long> A_inactive(n_inactive, n_dim_);
  nnz = n_bounds_ - n_xlb - n_xub;
  if (sub_.A_ub) nnz += sub_.A_ub->nonZeros() - A_active.nonZeros() + n_xlb + n_xub;
  A_inactive.reserve(nnz);
  n_xbound = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
    A_inactive.startVec(i);
    if (sub_.b_ub)
      for (Eigen::SparseMatrix<double>::InnerIterator it(*sub_.A_ub, i); it; ++it)
        if (!ub_activity.coeff(it.row()))
          A_inactive.insertBack(inactive_count_ub.coeff(it.row()), i) = it.value();
    if (sub_.lb && !xlb_activity[i] && sub_.lb->coeff(i) > -INF) {
      A_inactive.insertBack(n_ub_inactive + n_xbound, i) = -1;
      n_xbound++;
    }
    if (sub_.ub && !xub_activity[i] && sub_.ub->coeff(i) < INF) {
      A_inactive.insertBack(n_ub_inactive + n_xbound, i) = 1;
      n_xbound++;
    }
  }
  A_inactive.finalize();

  // Inactive b vector
  Eigen::VectorXd b_inactive(n_inactive);
  for (size_t i = 0; i < n_ub; ++i)
    if (!ub_activity.coeff(i))
      b_inactive.coeffRef(inactive_count_ub.coeff(i)) = sub_.b_ub->coeff(i);
  n_xbound = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
    if (sub_.lb && !xlb_activity[i] && sub_.lb->coeff(i) > -INF) {
      b_inactive.coeffRef(n_ub_inactive + n_xbound) = -sub_.lb->coeff(i);
      n_xbound++;
    }
    if (sub_.ub && !xub_activity[i] && sub_.ub->coeff(i) < INF) {
      b_inactive.coeffRef(n_ub_inactive + n_xbound) = sub_.ub->coeff(i);
      n_xbound++;
    }
  }

  // Inactive C matrix
  Eigen::SparseMatrix<double> C_inactive(n_inactive, n_dim_master_);
  nnz = 0;
  if (sub_.C_ub) nnz += sub_.C_ub->nonZeros();
  C_inactive.reserve(nnz);
  for (size_t i = 0; i < n_dim_master_; ++i) {
    C_inactive.startVec(i);
    if (sub_.b_ub)
      for (Eigen::SparseMatrix<double>::InnerIterator it(*sub_.C_ub, i); it; ++it)
        if (!ub_activity.coeff(it.row()))
          C_inactive.insertBack(inactive_count_ub.coeff(it.row()), i) = it.value();
  }
  C_inactive.finalize();

  // Determine critical region -----------------------------------------------------------------------------------------
  Eigen::VectorXd cc = sub_.C_ub->transpose() * res.dual_ub;
  CriticalRegion cr;
  cr.alpha = cc.transpose();
  cr.beta = res.fun - cc.dot(x);

  // Compute QR factorization
  // The inverse can be computed directly when the set of constraints are full rank. However, we ignore that case
  // because we need to determine rank and QR factorization is rank-revealing in Eigen. It allows us to decompose the
  // matrix which will result in faster inversion.
  Eigen::SPQR<Eigen::SparseMatrix<double>> active_solver(A_active);
  long n_rank = active_solver.rank();

  // Handle dual degeneracy
  if (n_rank < n_dim_) {
    cr.n_add = n_dim_ - n_rank;
    Eigen::SparseMatrix<double> U1 = active_solver.matrixR().topLeftCorner(n_rank, n_rank);
    Eigen::SparseMatrix<double> U2 = active_solver.matrixR().topRightCorner(n_rank, cr.n_add);

    Eigen::MatrixXd PD = active_solver.matrixQ().transpose() * C_active.toDense();
    Eigen::MatrixXd P = PD.topLeftCorner(n_rank, n_dim_master_);

    Eigen::VectorXd qr = active_solver.matrixQ().transpose() * b_active;
    Eigen::VectorXd q = qr.head(n_rank);

    Eigen::SparseMatrix<double> EF = A_inactive * active_solver.colsPermutation();
    Eigen::SparseMatrix<double> E = EF.topLeftCorner(n_inactive, n_rank);
    Eigen::SparseMatrix<double> F = EF.topRightCorner(n_inactive, cr.n_add);

    Eigen::TriangularView<const Eigen::SparseMatrix<double>, Eigen::Upper> U1_tri = U1.triangularView<Eigen::Upper>();
    Eigen::MatrixXd U1P = U1_tri.solve(P);
    Eigen::MatrixXd U1U2 = U1_tri.solve(U2.toDense());
    Eigen::VectorXd U1q = U1_tri.solve(q);

    cr.A = C_inactive - E * U1P;
    cr.b = b_inactive - E * U1q;

//    cr.c_add = (c_->transpose() * active_solver.colsPermutation()).tail(cr.n_add);
    cr.A_add = F - E * U1U2;

    cr.A_add.prune(1e-8);
  }
  else {
    Eigen::SparseMatrix<double> U1 = active_solver.matrixR().topLeftCorner(n_rank, n_rank);

    Eigen::VectorXd Qtb_active = (active_solver.matrixQ().transpose() * b_active).head(n_rank);
    Eigen::MatrixXd QtC_active = (active_solver.matrixQ().transpose() * C_active.toDense()).topRows(n_rank);

    Eigen::TriangularView<const Eigen::SparseMatrix<double>, Eigen::Upper> U1_tri = U1.triangularView<Eigen::Upper>();
    Eigen::MatrixXd U1QtC = U1_tri.solve(QtC_active);
    Eigen::VectorXd U1Qtb = U1_tri.solve(Qtb_active);

    Eigen::SparseMatrix<double> A_inactive_perm = A_inactive * active_solver.colsPermutation();
    cr.A = C_inactive - A_inactive_perm * U1QtC;
    cr.b = b_inactive - A_inactive_perm * U1Qtb;
    cr.n_add = 0;
  }

  cr.A.prune(1e-8);

  return cr;
}
