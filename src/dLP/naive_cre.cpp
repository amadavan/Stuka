//
// Created by Avinash on 8/31/2020.
//

#include <stuka/dLP/naive_cre.h>

stuka::dLP::NaiveCRE::NaiveCRE(const stuka::dLP::DecomposedLinearProgram &dlp,
                               const stuka::Options &opts) : BaseDLPSolver(dlp, opts), opts_(opts) {

  // Set properties
  first_run_ = true;
  n_sub_calls_ = 0;
  n_sub_ = dlp.c.size() - 1;
  n_dim_master_ = dlp.c.back()->size();
  opt_val_ = INF;
  n_active_ = 0;
  n_visit_ = 0;

  // Generate master problem
  master_lp_.c = dlp.c.back();
  master_lp_.A_ub = dlp.A_ub.back();
  master_lp_.b_ub = dlp.b_ub.back();
  master_lp_.A_eq = dlp.A_eq.back();
  master_lp_.b_eq = dlp.b_eq.back();
  master_lp_.lb = dlp.lb.back();
  master_lp_.ub = dlp.ub.back();

  // Set subproblems
  subproblems_.reserve(n_sub_);
  for (size_t i = 0; i < n_sub_; ++i) {
    Subproblem sub{dlp.c[i], dlp.A_ub[i], dlp.b_ub[i], dlp.C_ub[i], dlp.A_eq[i], dlp.b_eq[i], dlp.C_eq[i], dlp.lb[i],
                   dlp.ub[i]};

    subproblems_.push_back({NaiveCRESubproblem(std::move(sub), opts_)});
  }

  // Set useful constants
  n_eq_ = (master_lp_.b_eq) ? master_lp_.b_eq->size() : 0;
  n_ub_ = (master_lp_.b_ub) ? master_lp_.b_ub->size() : 0;
  n_ub_cr_ = 0;

  bone_eq_ = Eigen::Matrix<bool, Eigen::Dynamic, 1>::Constant(n_eq_, true);

  eye_ = Eigen::SparseMatrix<double>(n_dim_master_, n_dim_master_);
  eye_.setIdentity();
  negative_eye_ = -eye_;

  // Set current cost to cost of master problem
  cr_.alpha = 0 * *master_lp_.c;
  cr_.beta = 0;
  cr_.n_add = 0;
  critical_regions_ = std::vector<CriticalRegionData>(n_sub_);

  // Set initial point
  x_ = Eigen::VectorXd::Constant(n_dim_master_, 1);
  if (opts_.x0.size() == n_dim_master_) {
    x_ = opts_.x0;
  } else if (master_lp_.c->maxCoeff() > 1e-8 && master_lp_.c->minCoeff() < -1e-8) {
    OptimizeState res = util::createSolver(master_lp_, opts_)->solve();
    if (res.status == 2) x_ = res.x;
  }

  // Create residual solver
  residual_qp_.Q = std::make_shared<Eigen::SparseMatrix<double>>(n_dim_master_, n_dim_master_);
  residual_qp_.Q->setIdentity();
  residual_qp_.A_eq = std::make_shared<Eigen::SparseMatrix<double>>(n_dim_master_, n_dim_master_);
  residual_qp_.A_eq->setIdentity();
  (*residual_qp_.A_eq) *= -1;
  residual_qp_.A_eq->conservativeResize(n_dim_master_ + 1, n_dim_master_);
  residual_qp_.b_eq = std::make_shared<Eigen::VectorXd>(n_dim_master_ + 1);
  residual_qp_.b_eq->setZero();
  residual_qp_.b_eq->coeffRef(n_dim_master_) = 1.;

}

void stuka::dLP::NaiveCRE::iterate() {

  for (size_t i = 0; i < n_sub_; ++i) {
    if (!critical_regions_[i].cr.in(x_)) {
      n_sub_calls_++;

      // Remove current critical region from problem
      if (!first_run_) removeCriticalRegion(critical_regions_[i]);

      critical_regions_[i].cr = subproblems_[i].computeCriticalRegion(x_);

      addCriticalRegion(critical_regions_[i]);
    }
  }

  first_run_ = false;

  OptimizeState res = solveMasterProblem();
  res.fun += cr_.beta;

  // Remove any additional variables
  Eigen::VectorXd x = res.x.head(n_dim_master_);

  if (res.fun + 1e-8 < opt_val_) {
    x_opt_ = x;
    opt_val_ = res.fun;
    a_opt_ = computeResidualFirst(res);
  } else {
    a_opt_ = computeResidualSubsequent();
  }

  // Set next iterate
  x_ = x - a_opt_ / a_opt_.norm() * opts_.cre_step;

}

bool stuka::dLP::NaiveCRE::terminate() {
  return (a_opt_.squaredNorm() < opts_.tol);
}

const stuka::OptimizeState stuka::dLP::NaiveCRE::getState() {
  OptimizeState res;

  res.fun = (*master_lp_.c + cr_.alpha).transpose() * x_opt_ + cr_.beta;
  res.x = x_opt_;
  res.error = a_opt_.squaredNorm();
  res.status = (res.error < opts_.tol) ? 2 : 1;
  res.nit_sub = n_sub_calls_;

  return res;
}

void stuka::dLP::NaiveCRE::removeCriticalRegion(const stuka::dLP::NaiveCRE::CriticalRegionData &crdata) {

  // Useful variables
  const CriticalRegion &cr = crdata.cr;
  long n_con_cr = cr.b.size();
  long n_con_cr_ = cr_.b.size();

  mtx_cre_.lock(); // Need a mutex when updating global properties

  // Remove critical region from aggregate cost
  cr_.alpha -= cr.alpha;
  cr_.beta -= cr.beta;

  // Remove constraints
  Eigen::SparseMatrix<double> A = Eigen::SparseMatrix<double>(n_con_cr_ - n_con_cr, n_dim_master_);
  A.reserve(cr_.A.nonZeros());
  for (size_t i = 0; i < n_dim_master_; ++i) {
    A.startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(cr_.A, i); it; ++it) {
      if (it.row() < crdata.ind_constr)
        A.insertBack(it.row(), i) = it.value();
      else if (it.row() > crdata.ind_constr + n_con_cr)
        A.insertBack(it.row() - n_con_cr, i) = it.value();
    }
  }
  A.finalize();
  cr_.A = A;

  // Remove additional variables
  Eigen::SparseMatrix<double> A_add = Eigen::SparseMatrix<double>(n_con_cr_ - n_con_cr, cr_.n_add - cr.n_add);
  A_add.reserve(cr_.A_add.nonZeros());
  for (size_t i = 0; i < crdata.ind_add; ++i) {
    A_add.startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(cr_.A_add, i); it; ++it) {
      if (it.row() < crdata.ind_constr)
        A_add.insertBack(it.row(), i) = it.value();
      else if (it.row() > crdata.ind_constr + n_con_cr)
        A_add.insertBack(it.row() - n_con_cr, i) = it.value();
    }
  }
  for (size_t i = crdata.ind_add + cr.n_add; i < cr_.n_add; ++i) {
    A_add.startVec(i - cr.n_add);
    for (Eigen::SparseMatrix<double>::InnerIterator it(cr_.A_add, i); it; ++it) {
      if (it.row() < crdata.ind_constr)
        A_add.insertBack(it.row(), i - cr.n_add) = it.value();
      else if (it.row() > crdata.ind_constr + n_con_cr)
        A_add.insertBack(it.row() - n_con_cr, i - cr.n_add) = it.value();

    }
  }
  A_add.finalize();
  cr_.A_add = A_add;

  // Remove constraint bounds
  Eigen::VectorXd b = Eigen::VectorXd(n_con_cr_ - n_con_cr);
  b.head(crdata.ind_constr) = cr_.b.head(crdata.ind_constr);
  b.tail(n_con_cr_ - n_con_cr - crdata.ind_constr) = cr_.b.tail(n_con_cr_ - n_con_cr - crdata.ind_constr);
  cr_.b = b;

  n_ub_cr_ -= n_con_cr;
  cr_.n_add -= cr.n_add;

  // Shift indices of all other critical regions if they appear after
  for (size_t j = 0; j < n_sub_; ++j) {
    CriticalRegionData &crdata2 = critical_regions_[j];
    if (&crdata2 == &crdata) continue;
    if (crdata.ind_constr < crdata2.ind_constr)
      crdata2.ind_constr -= n_con_cr;
    if (crdata.ind_add < crdata2.ind_add)
      crdata2.ind_add -= cr.n_add;
  }

  mtx_cre_.unlock();
}

void stuka::dLP::NaiveCRE::addCriticalRegion(stuka::dLP::NaiveCRE::CriticalRegionData &crdata) {

  const CriticalRegion &cr = crdata.cr;

  mtx_cre_.lock(); // Need a mutex when updating global properties

  const size_t n_con_cr_ = cr_.A.rows();

  // Add cost of critical region to aggregate cost
  cr_.alpha += cr.alpha;
  cr_.beta += cr.beta;

  // Add constraints associated with critical region
  Eigen::SparseMatrix<double> A = Eigen::SparseMatrix<double>(n_con_cr_ + cr.A.rows(), n_dim_master_);
  A.reserve(cr_.A.nonZeros() + cr.A.nonZeros());
  // nominal variables
  for (size_t i = 0; i < n_dim_master_; ++i) {
    A.startVec(i);
    if (n_con_cr_ > 0)
      for (Eigen::SparseMatrix<double>::InnerIterator it(cr_.A, i); it; ++it)
        A.insertBack(it.row(), i) = it.value();
    for (Eigen::SparseMatrix<double>::InnerIterator it(cr.A, i); it; ++it)
      A.insertBack(n_con_cr_ + it.row(), i) = it.value();
  }
  A.finalize();

  // Add additional variables
  Eigen::SparseMatrix<double> A_add = Eigen::SparseMatrix<double>(n_con_cr_ + cr.A.rows(), cr_.n_add + cr.n_add);
  A_add.reserve(cr_.A_add.nonZeros() + cr.A_add.nonZeros());
  // existing additional variables
  for (size_t i = 0; i < cr_.n_add; ++i) {
    A_add.startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(cr_.A_add, i); it; ++it)
      A_add.insertBack(it.row(), i) = it.value();
  }
  // new additional variables
  for (size_t i = 0; i < cr.n_add; ++i) {
    A_add.startVec(cr_.n_add + i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(cr.A_add, i); it; ++it)
      A_add.insertBack(n_con_cr_ + it.row(), cr_.n_add + i) = it.value();
  }
  A_add.finalize();

  cr_.A = A;
  cr_.A_add = A_add;

  // Constraint bounds
  cr_.b.conservativeResize(cr_.b.size() + cr.b.size());
  cr_.b.tail(cr.b.size()) = cr.b;

  // Add new additional variables if necessary
  if (cr.n_add > 0) {
    // Set additional variables index and update total constraints
    crdata.ind_add = cr_.n_add;
    cr_.n_add += cr.n_add;
  }

  // Set constraint index and update total constraints
  crdata.ind_constr = n_ub_ + n_ub_cr_;
  n_ub_cr_ += cr.b.size();

  mtx_cre_.unlock();
}

stuka::OptimizeState stuka::dLP::NaiveCRE::solveMasterProblem() {

  LP::LinearProgram prog(master_lp_);

  if (cr_.n_add > 0) {
    std::shared_ptr<Eigen::VectorXd> c = std::make_shared<Eigen::VectorXd>(*master_lp_.c + cr_.alpha);
    c->conservativeResize(n_dim_master_ + cr_.n_add);
    for (size_t i = 0; i < cr_.n_add; ++i) c->coeffRef(n_dim_master_ + i) = 0.;
    prog.c = c;

    if (master_lp_.A_eq) {
      std::shared_ptr<Eigen::SparseMatrix<double>>
          A_eq = std::make_shared<Eigen::SparseMatrix<double>>(*master_lp_.A_eq);
      A_eq->conservativeResize(A_eq->rows(), A_eq->cols() + cr_.n_add);
      prog.A_eq = A_eq;
    }

    if (master_lp_.lb) {
      std::shared_ptr<Eigen::VectorXd> lb = std::make_shared<Eigen::VectorXd>(*master_lp_.lb);
      lb->conservativeResize(n_dim_master_ + cr_.n_add);
      for (size_t i = 0; i < cr_.n_add; ++i) lb->coeffRef(n_dim_master_ + i) = -INF;
      prog.lb = lb;
    }

    if (master_lp_.ub) {
      std::shared_ptr<Eigen::VectorXd> ub = std::make_shared<Eigen::VectorXd>(*master_lp_.ub);
      ub->conservativeResize(n_dim_master_ + cr_.n_add);
      for (size_t i = 0; i < cr_.n_add; ++i) ub->coeffRef(n_dim_master_ + i) = INF;
      prog.ub = ub;
    }
  }

  // Combine critical region and master problem
  std::shared_ptr<Eigen::SparseMatrix<double>>
      A = std::make_shared<Eigen::SparseMatrix<double>>(n_ub_ + cr_.A.rows(), n_dim_master_ + cr_.n_add);
  if (master_lp_.A_ub) A->reserve(master_lp_.A_ub->nonZeros() + cr_.A.nonZeros() + cr_.A_add.nonZeros());
  else A->reserve(cr_.A.nonZeros() + cr_.A.nonZeros() + cr_.A_add.nonZeros());
  // nominal variables
  for (size_t i = 0; i < n_dim_master_; ++i) {
    A->startVec(i);
    if (master_lp_.A_ub)
      for (Eigen::SparseMatrix<double>::InnerIterator it(*master_lp_.A_ub, i); it; ++it)
        A->insertBack(it.row(), i) = it.value();
    for (Eigen::SparseMatrix<double>::InnerIterator it(cr_.A, i); it; ++it)
      A->insertBack(n_ub_ + it.row(), i) = it.value();
  }
  // additional variables
  for (size_t i = 0; i < cr_.n_add; ++i) {
    A->startVec(n_dim_master_ + i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(cr_.A_add, i); it; ++it)
      A->insertBack(n_ub_ + it.row(), n_dim_master_ + i) = it.value();
  }
  A->finalize();
  prog.A_ub = A;

  std::shared_ptr<Eigen::VectorXd> b = std::make_shared<Eigen::VectorXd>(n_ub_ + cr_.A.rows());
  if (master_lp_.b_ub) b->head(n_ub_) = *master_lp_.b_ub;
  b->tail(cr_.A.rows()) = cr_.b;
  prog.b_ub = b;

  std::unique_ptr<LP::BaseLPSolver> solver = util::createSolver(prog, opts_);

  OptimizeState res;
  try {
    res = solver->solve();
  } catch (std::exception &e) {
//    std::cout << "Failed to solve master problem: " << res.status << std::endl;
//    std::cout << e.getErrorCode() << "\t" << e.getMessage() << std::endl;
    throw std::runtime_error("iterate: unable to solve master problem");
  }
  if (res.status != 2) {
    std::cout << "Failed to solve master problem: " << res.status << std::endl;
    throw std::runtime_error("iterate: unable to solve master problem");
  }

  return res;
}

// TODO: Intelligently remove constant variables from solution of master problem
const Eigen::VectorXd
stuka::dLP::NaiveCRE::computeResidualFirst(const OptimizeState &res) {

  // Remove additional variables from residual problem
//  if (n_active_ + n_visit_ > 0) residual_solver_->getQP().removeBackVars(n_active_ + n_visit_);
  n_visit_ = 1;

  // Determine active constraints
  util::DenseOps::ActiveSet ub_activity = util::DenseOps::ActiveSet::Constant(n_ub_, false);
  util::DenseOps::ActiveSet xlb_activity = util::DenseOps::ActiveSet::Constant(n_dim_master_, false);
  util::DenseOps::ActiveSet xub_activity = util::DenseOps::ActiveSet::Constant(n_dim_master_, false);

  // Construct set of active constraints
  std::vector<std::pair<const Eigen::SparseMatrix<double>, const util::DenseOps::ActiveSet>> A_active_l =
      std::vector<std::pair<const Eigen::SparseMatrix<double>, const util::DenseOps::ActiveSet>>();

  if (master_lp_.A_eq) A_active_l.emplace_back(*master_lp_.A_eq, bone_eq_);
  if (master_lp_.b_ub) {
    ub_activity = res.dual_ub.head(n_ub_).array() > 0;
    A_active_l.emplace_back(*master_lp_.A_ub, ub_activity);
  }
  if (master_lp_.lb) {
    xlb_activity = res.dual_x_lb.head(n_dim_master_).array() > 0;
    A_active_l.emplace_back(negative_eye_, xlb_activity);
  }
  if (master_lp_.ub) {
    xub_activity = res.dual_x_ub.head(n_dim_master_).array() > 0;
    A_active_l.emplace_back(eye_, xub_activity);
  }

  // Create useful constants
  long n_ub = ub_activity.count();
  long n_xlb = xlb_activity.count();
  long n_xub = xub_activity.count();
  n_active_ = n_eq_ + n_ub + n_xlb + n_xub;

  Eigen::SparseMatrix<double> A_active = util::SparseOps::vstack_rows<double>(A_active_l);

  A_activeT_ = A_active.transpose();
  previous_directions_ = Eigen::MatrixXd(n_dim_master_, 1);
  previous_directions_.col(0) = *master_lp_.c + cr_.alpha;

  // Determine gradient adjustment
  Eigen::SparseMatrix<double> grad_adjust = Eigen::SparseMatrix<double>(n_dim_master_, n_dim_master_);
  grad_adjust.setIdentity();
  if (n_active_ > 0) {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> AAt;
    AAt.compute(A_active * A_active.transpose());
    grad_adjust -= A_active.transpose() * AAt.solve(A_active);
    grad_adjust.prune(1e-8);
  }

  return grad_adjust * (*master_lp_.c + cr_.alpha);
}

const Eigen::VectorXd stuka::dLP::NaiveCRE::computeResidualSubsequent() {
  n_visit_++;

  previous_directions_.conservativeResize(n_dim_master_, n_visit_);
  previous_directions_.col(n_visit_ - 1) = *master_lp_.c + cr_.alpha;

  QP::QuadraticProgram prog;
  prog.Q = std::make_shared<Eigen::SparseMatrix<double>>(*residual_qp_.Q);
  prog.Q->conservativeResize(n_dim_master_ + n_active_ + n_visit_, n_dim_master_ + n_active_ + n_visit_);

  prog.A_eq = std::make_shared<Eigen::SparseMatrix<double>>(n_dim_master_ + 1, n_dim_master_ + n_active_ + n_visit_);
  prog.A_eq->conservativeResize(n_dim_master_ + 1, n_dim_master_ + n_active_ + n_visit_);
  prog.A_eq->reserve(n_dim_master_ + A_activeT_.nonZeros() + previous_directions_.nonZeros() + n_active_);
  for (size_t i = 0; i < n_dim_master_; ++i) {
    prog.A_eq->startVec(i);
    prog.A_eq->insertBack(i, i) = -1;
  }
  for (size_t i = 0; i < n_active_; ++i) {
    prog.A_eq->startVec(n_dim_master_ + i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(A_activeT_, i); it; ++it)
      prog.A_eq->insertBack(it.row(), n_dim_master_ + i) = it.value();
  }
  for (size_t i = 0; i < n_visit_; ++i) {
    prog.A_eq->startVec(n_dim_master_ + n_active_ + i);
    for (Eigen::MatrixXd::InnerIterator it(previous_directions_, i); it; ++it) {
      if (abs(it.value()) > 1e-6)
        prog.A_eq->insertBack(it.row(), n_dim_master_ + n_active_ + i) = it.value();
    }
    prog.A_eq->insertBack(n_dim_master_, n_dim_master_ + n_active_ + i) = 1;
  }
  prog.A_eq->finalize();

  prog.b_eq = std::make_shared<Eigen::VectorXd>(n_dim_master_ + 1);
  prog.b_eq->setZero();
  prog.b_eq->coeffRef(n_dim_master_) = 1.;

  prog.lb = std::make_shared<Eigen::VectorXd>(n_dim_master_ + n_active_ + n_visit_);
  prog.lb->setConstant(-INF);
  for (size_t i = 0; i < n_active_ + n_visit_; ++i) prog.lb->coeffRef(n_dim_master_ + i) = 0;

  std::unique_ptr<QP::BaseQPSolver> residual_solver = util::createSolver(prog, opts_);

  OptimizeState res;
  try {
    res = residual_solver->solve();
  } catch (std::exception &e) {
//    std::cout << "Failed to solve residual problem: " << res.status << std::endl;
//    std::cout << e.getErrorCode() << "\t" << e.getMessage() << std::endl;
    throw std::runtime_error("iterate: unable to solve residual problem");
  }
  if (res.status != 2) {
    std::cout << "Failed to solve residual problem: " << res.status << std::endl;
    throw std::runtime_error("iterate: unable to solve residual problem");
  }

  return res.x.head(n_dim_master_);
}
