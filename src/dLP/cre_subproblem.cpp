//
// Created by Avinash Madavan on 1/15/19.
//

#include <stuka/dLP/cre_subproblem.h>

stuka::dLP::CRESubproblem::CRESubproblem(const stuka::dLP::Subproblem sub, const stuka::Options opts) : sub_(sub) {
  n_dim_ = sub_.c->size();
  n_dim_master_ = (sub_.C_ub) ? sub_.C_ub->cols() : (sub_.C_eq) ? sub_.C_eq->cols() : 0;

  bone_eq_ = util::DenseOps::ActiveSet::Constant((sub_.b_eq) ? sub_.b_eq->size() : 0, true);

  n_bounds_ = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
    if (sub_.lb && sub_.lb->coeff(i) > -INF) n_bounds_++;
    if (sub_.ub && sub_.ub->coeff(i) < INF) n_bounds_++;
  }

  b_lb_ = Eigen::VectorXd::Constant(n_dim_, INF);
  b_ub_ = Eigen::VectorXd::Constant(n_dim_, INF);

  xlb_constraints_ = util::DenseOps::ActiveSet::Constant(n_dim_, false);
  xub_constraints_ = util::DenseOps::ActiveSet::Constant(n_dim_, false);

  if (sub_.lb) {
    A_xlb_ = Eigen::SparseMatrix<double>(n_dim_, n_dim_);
    A_xlb_.setIdentity();
    A_xlb_ *= -1;
    b_lb_ = *sub_.lb;
    b_lb_ *= -1;

    for (size_t i = 0; i < n_dim_; ++i)
      xlb_constraints_.coeffRef(i) = sub_.lb->coeff(i) > -INF;

    A_xlb_ = util::SparseOps::get_rows(A_xlb_, xlb_constraints_);
    b_lb_ = util::DenseOps::get_rows(b_lb_, xlb_constraints_);
  }

  if (sub_.ub) {
    A_xub_ = Eigen::SparseMatrix<double>(n_dim_, n_dim_);
    A_xub_.setIdentity();
    b_ub_ = *sub_.ub;

    for (size_t i = 0; i < n_dim_; ++i)
      xub_constraints_.coeffRef(i) = sub_.ub->coeff(i) < INF;

    A_xub_ = util::SparseOps::get_rows(A_xub_, xub_constraints_);
    b_ub_ = util::DenseOps::get_rows(b_ub_, xub_constraints_);
  }

  LP::LinearProgram lp;
  lp.c = sub_.c;
  lp.A_ub = sub_.A_ub;
  lp.b_ub = sub_.b_ub;
  lp.A_eq = sub_.A_eq;
  lp.b_eq = sub_.b_eq;
  lp.lb = sub_.lb;
  lp.ub = sub_.ub;

  solver_ = util::createSolver(lp, opts);
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

  // Construct list of matrices with subset of rows to include
  std::vector<util::SparseOps::MatRowPairD> A_active_l = std::vector<util::SparseOps::MatRowPairD>();
  std::vector<util::DenseOps::VecRowPairD> b_active_l = std::vector<util::DenseOps::VecRowPairD>();
  std::vector<util::SparseOps::MatRowPairD> C_active_l = std::vector<util::SparseOps::MatRowPairD>();
  std::vector<util::SparseOps::MatRowPairD> A_inactive_l = std::vector<util::SparseOps::MatRowPairD>();
  std::vector<util::DenseOps::VecRowPairD> b_inactive_l = std::vector<util::DenseOps::VecRowPairD>();
  std::vector<util::SparseOps::MatRowPairD> C_inactive_l = std::vector<util::SparseOps::MatRowPairD>();

  if (sub_.A_eq) {
    A_active_l.push_back(std::make_pair(*sub_.A_eq, bone_eq_));
    b_active_l.push_back(std::make_pair(*sub_.b_eq, bone_eq_));
    C_active_l.push_back(std::make_pair(*sub_.C_eq, bone_eq_));
  }
  if (sub_.A_ub) {
    util::DenseOps::ActiveSet ub_activity = res.dual_ub.array() > 0;
    util::DenseOps::ActiveSet ub_inactivity = res.dual_ub.array() <= 0;

    A_active_l.push_back(std::make_pair(*sub_.A_ub, ub_activity));
    b_active_l.push_back(std::make_pair(*sub_.b_ub, ub_activity));
    C_active_l.push_back(std::make_pair(*sub_.C_ub, ub_activity));

    A_inactive_l.push_back(std::make_pair(*sub_.A_ub, ub_inactivity));
    b_inactive_l.push_back(std::make_pair(*sub_.b_ub, ub_inactivity));
    C_inactive_l.push_back(std::make_pair(*sub_.C_ub, ub_inactivity));
  }
  if (sub_.lb) {
    util::DenseOps::ActiveSet xlb_activity = res.dual_x_lb.array() > 0;
    util::DenseOps::ActiveSet xlb_inactivity = res.dual_x_lb.array() <= 0;

    xlb_activity = util::DenseOps::get_rows<bool, Eigen::Dynamic, 1>(xlb_activity, xlb_constraints_);
    xlb_inactivity = util::DenseOps::get_rows<bool, Eigen::Dynamic, 1>(xlb_inactivity, xlb_constraints_);

    if (xlb_activity.size() > 0 && xlb_activity.count() > 0) {
      A_active_l.push_back(std::make_pair(A_xlb_, xlb_activity));
      b_active_l.push_back(std::make_pair(b_lb_, xlb_activity));
    }

    if (xlb_inactivity.size() > 0 && xlb_inactivity.count() > 0) {
      A_inactive_l.push_back(std::make_pair(A_xlb_, xlb_inactivity));
      b_inactive_l.push_back(std::make_pair(b_lb_, xlb_inactivity));
    }
  }
  if (sub_.ub) {
    util::DenseOps::ActiveSet xub_activity = res.dual_x_ub.array() > 0;
    util::DenseOps::ActiveSet xub_inactivity = res.dual_x_ub.array() <= 0;

    xub_activity = util::DenseOps::get_rows<bool, Eigen::Dynamic, 1>(xub_activity, xub_constraints_);
    xub_inactivity = util::DenseOps::get_rows<bool, Eigen::Dynamic, 1>(xub_inactivity, xub_constraints_);

    if (xub_activity.size() > 0 && xub_activity.count() > 0) {
      A_active_l.push_back(std::make_pair(A_xub_, xub_activity));
      b_active_l.push_back(std::make_pair(b_ub_, xub_activity));
    }

    if (xub_inactivity.size() > 0 && xub_inactivity.count() > 0) {
      A_inactive_l.push_back(std::make_pair(A_xub_, xub_inactivity));
      b_inactive_l.push_back(std::make_pair(b_ub_, xub_inactivity));
    }
  }

  // Construct set of active constraints -------------------------------------------------------------------------------
  Eigen::SparseMatrix<double> A_active = util::SparseOps::vstack_rows<double>(A_active_l);        // Active A matrix
  Eigen::VectorXd b_active = util::DenseOps::vstack_rows<double, Eigen::Dynamic, 1>(b_active_l);  // Active b vector
  Eigen::SparseMatrix<double> C_active = util::SparseOps::vstack_rows<double>(C_active_l);        // Active C matrix
  C_active.conservativeResize(A_active.rows(), n_dim_master_);

  // Construct set of inactive constraints -----------------------------------------------------------------------------
  Eigen::SparseMatrix<double, 0, long> A_inactive = util::SparseOps::vstack_rows(A_inactive_l);   // Inactive A matrix
  Eigen::VectorXd b_inactive = util::DenseOps::vstack_rows(b_inactive_l);                         // Inactive b vector

  // Inactive C matrix
  Eigen::SparseMatrix<double> C_inactive = util::SparseOps::vstack_rows<double>(C_inactive_l);
  C_inactive.conservativeResize(A_inactive.rows(), n_dim_master_);

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
  cr.n_add = n_dim_ - n_rank;

  // Handle dual degeneracy
  if (n_rank < n_dim_) {
    Eigen::SparseMatrix<double> U1 = active_solver.matrixR().topLeftCorner(n_rank, n_rank);
    Eigen::SparseMatrix<double> U2 = active_solver.matrixR().topRightCorner(n_rank, cr.n_add);

    Eigen::MatrixXd PD = active_solver.matrixQ().transpose() * C_active.toDense();
    Eigen::MatrixXd P = PD.topLeftCorner(n_rank, n_dim_master_);

    Eigen::VectorXd qr = active_solver.matrixQ().transpose() * b_active;
    Eigen::VectorXd q = qr.head(n_rank);

    size_t n_inactive = A_inactive.rows();
    Eigen::SparseMatrix<double> EF = A_inactive * active_solver.colsPermutation();
    Eigen::SparseMatrix<double> E = EF.topLeftCorner(n_inactive, n_rank);
    Eigen::SparseMatrix<double> F = EF.topRightCorner(n_inactive, cr.n_add);

    Eigen::TriangularView<const Eigen::SparseMatrix<double>, Eigen::Upper> U1_tri = U1.triangularView<Eigen::Upper>();
    Eigen::MatrixXd U1P = U1_tri.solve(P);
    Eigen::MatrixXd U1U2 = U1_tri.solve(U2.toDense());
    Eigen::VectorXd U1q = U1_tri.solve(q);

    cr.A = C_inactive - E * U1P;
    cr.b = b_inactive - E * U1q;

    cr.A_add = F - E * U1U2;

    cr.A_add.prune(1e-8);
  } else {
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
