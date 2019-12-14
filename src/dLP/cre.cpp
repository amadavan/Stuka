//
// Created by Avinash Madavan on 1/21/19.
//

#include <stuka/dLP/cre.h>

stuka::dLP::CRE::CRE(const stuka::dLP::DecomposedLinearProgram &dlp,
                     const stuka::Options &opts) : BaseDLPSolver(dlp, opts) {

  // Set properties
  first_run_ = true;
  n_sub_calls_ = 0;
  n_sub_ = dlp.c.size() - 1;
  n_dim_master_ = dlp.c.back()->size();
  opt_val_ = INF;
  n_active_ = 0;
  n_visit_ = 0;
  n_add_ = 0;
  n_con_ub_ = (master_lp_.b_ub) ? master_lp_.b_ub->size() : 0;

  // Generate master problem
  master_lp_.c = dlp.c.back();
  master_lp_.A_ub = dlp.A_ub.back();
  master_lp_.b_ub = dlp.b_ub.back();
  master_lp_.A_eq = dlp.A_eq.back();
  master_lp_.b_eq = dlp.b_eq.back();
  master_lp_.lb = dlp.lb.back();
  master_lp_.ub = dlp.ub.back();

  master_solver_ = util::createSolver(master_lp_, opts);

  // Set subproblems
  subproblems_.reserve(n_sub_);
  for (size_t i = 0; i < n_sub_; ++i) {
    Subproblem sub{dlp.c[i], dlp.A_ub[i], dlp.b_ub[i], dlp.C_ub[i], dlp.A_eq[i], dlp.b_eq[i], dlp.C_eq[i], dlp.lb[i],
                   dlp.ub[i]};

    subproblems_.push_back({CRESubproblem(std::move(sub))});
  }

  // Set current cost to cost of master problem
  cost_alpha_ = *master_lp_.c;
  cost_beta_ = 0;
  critical_regions_ = std::vector<CriticalRegionData>(n_sub_);

  // Set initial point
  x_ = Eigen::VectorXd(n_dim_master_);
  x_.setConstant(1.);
  if (opts.x0.size() == n_dim_master_) {
    x_ = opts.x0;
  } else if (master_lp_.c->maxCoeff() > 1e-8 && master_lp_.c->minCoeff() < -1e-8) {
    OptimizeState res = util::createSolver(master_lp_, opts)->solve();
    if (res.status == 2) x_ = res.x;
  }

  // Create residual solver
  QP::QuadraticProgram residual;
  residual.Q = std::make_shared<Eigen::SparseMatrix<double>>(n_dim_master_, n_dim_master_);
  residual.Q->setIdentity();
  residual.A_eq = std::make_shared<Eigen::SparseMatrix<double>>(n_dim_master_, n_dim_master_);
  residual.A_eq->setIdentity();
  (*residual.A_eq) *= -1;
  residual.A_eq->conservativeResize(n_dim_master_ + 1, n_dim_master_);
  residual.b_eq = std::make_shared<Eigen::VectorXd>(n_dim_master_ + 1);
  residual.b_eq->setZero();
  residual.b_eq->coeffRef(n_dim_master_) = 1.;
  residual_solver_ = util::createSolver(residual, opts);

}

void stuka::dLP::CRE::iterate() {

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

  // Remove any additional variables
  Eigen::VectorXd x = res.x.head(n_dim_master_);

  if (res.fun + cost_beta_ + 1e-8 < opt_val_) {
    x_opt_ = x;
    opt_val_ = res.fun + cost_beta_;
    a_opt_ = computeResidualFirst(res);
  } else {
    a_opt_ = computeResidualSubsequent();
  }

  // Set next iterate
  x_ = x - a_opt_ / a_opt_.norm() * 1e-4;

}

bool stuka::dLP::CRE::terminate() {
  // TODO: set tolerance from options
  return (a_opt_.squaredNorm() < 1e-8);
}

const stuka::OptimizeState stuka::dLP::CRE::getState() {
  OptimizeState res;

  res.x = x_opt_;
  res.error = a_opt_.squaredNorm();
  res.status = (res.error < 1e-8) ? 2 : 1;
  res.nit_sub = n_sub_calls_;

  return res;
}

void stuka::dLP::CRE::removeCriticalRegion(const stuka::dLP::CRE::CriticalRegionData &crdata) {
  mtx_cre_.lock(); // Need a mutex when updating global properties

  // Useful variables
  const CriticalRegion &cr = crdata.cr;
  long n_con_cr_ = cr.b.size();

  // Remove critical region from aggregate cost
  cost_alpha_ -= cr.alpha;
  cost_beta_ -= cr.beta;

  master_solver_->getLP().removeConstrs_ub(crdata.ind_constr, n_con_cr_);
  n_con_ub_ -= n_con_cr_;

  // Remove additional variables and constraints
  master_solver_->getLP().removeVars(n_dim_master_ + crdata.ind_add, cr.n_add);
  n_add_ -= cr.n_add;

  // Shift indices of all other critical regions if they appear after
  for (size_t j = 0; j < n_sub_; ++j) {
    CriticalRegionData &crdata2 = critical_regions_[j];
    if (&crdata2 == &crdata) continue;
    if (crdata.ind_constr < crdata2.ind_constr)
      crdata2.ind_constr -= n_con_cr_;
    if (crdata.ind_add < crdata2.ind_add)
      crdata2.ind_add -= cr.n_add;
  }

  mtx_cre_.unlock();
}

void stuka::dLP::CRE::addCriticalRegion(stuka::dLP::CRE::CriticalRegionData &crdata) {
  mtx_cre_.lock(); // Need a mutex when updating global properties

  CriticalRegion &cr = crdata.cr;

  // Add cost of critical region to aggregate cost
  cost_alpha_ += cr.alpha;
  cost_beta_ += cr.beta;

  // Add constraints associated with critical region
  std::shared_ptr<Eigen::SparseMatrix<double>> A = std::make_shared<Eigen::SparseMatrix<double>>(cr.A);
  A->conservativeResize(A->rows(), n_dim_master_ + n_add_);
  master_solver_->getLP().addConstrs_ub(A, std::make_shared<Eigen::VectorXd>(cr.b));

  // Add new additional variables if necessary
  if (cr.n_add > 0) {
    std::shared_ptr<Eigen::SparseMatrix<double>> A_add =
        std::make_shared<Eigen::SparseMatrix<double>>(n_con_ub_ + cr.b.size(), cr.n_add);
    A_add->reserve(cr.A_add.nonZeros());
    for (size_t i = 0; i < cr.n_add; ++i) {
      A_add->startVec(i);
      for (Eigen::SparseMatrix<double>::InnerIterator itAdd(cr.A_add, i); itAdd; ++itAdd)
        A_add->insertBack(n_con_ub_ + itAdd.row(), i) = itAdd.value();
    }
    A_add->finalize();

    std::shared_ptr<Eigen::VectorXd> c_add = std::make_shared<Eigen::VectorXd>(cr.n_add);
    c_add->setZero();

    master_solver_->getLP().addVars(nullptr, A_add, nullptr, nullptr, nullptr);

    // Set additional variables index and update total constraints
    crdata.ind_add = n_add_;
    n_add_ += cr.n_add;
  }

  // Set constraint index and update total constraints
  crdata.ind_constr = n_con_ub_;
  n_con_ub_ += cr.b.size();

  mtx_cre_.unlock();
}

stuka::OptimizeState stuka::dLP::CRE::solveMasterProblem() {

  std::shared_ptr<Eigen::VectorXd> c = std::make_shared<Eigen::VectorXd>(cost_alpha_);
  c->conservativeResize(n_dim_master_ + n_add_);
  // TODO: check why we need to initialize these values - should be uninitialized by default
  for (size_t i = 0; i < n_add_; ++i) c->coeffRef(n_dim_master_ + i) = 0.;
  master_solver_->getLP().setObjective(c);

  OptimizeState res;
  try {
    res = master_solver_->solve();
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
stuka::dLP::CRE::computeResidualFirst(const OptimizeState &res) {

  // Remove additional variables from residual problem
  if (n_active_ + n_visit_ > 0) residual_solver_->getQP().removeBackVars(n_active_ + n_visit_);
  n_visit_ = 1;

  long n_con_ub = (master_lp_.b_ub) ? master_lp_.b_ub->size() : 0;

  // Determine active constraints
  Eigen::Matrix<bool, Eigen::Dynamic, 1> ub_activity(n_con_ub);
  if (master_lp_.b_ub) ub_activity = res.dual_ub.head(n_con_ub).array() > GUROBI_TOLERANCE;
//  Eigen::Matrix<bool, Eigen::Dynamic, 1> ub_activity = dual_ub.array() > 0;
  Eigen::Matrix<bool, Eigen::Dynamic, 1> xlb_activity(n_dim_master_);
  Eigen::Matrix<bool, Eigen::Dynamic, 1> xub_activity(n_dim_master_);
  xlb_activity.setConstant(false);
  xub_activity.setConstant(false);
  if (master_lp_.lb) xlb_activity = res.dual_x_lb.head(n_dim_master_).array() > GUROBI_TOLERANCE;
  if (master_lp_.ub) xub_activity = res.dual_x_ub.head(n_dim_master_).array() > GUROBI_TOLERANCE;
//  for (size_t i = 0; i < n_dim_master_; ++i) {
//    if (master_lp_.lb && master_lp_.lb->coeff(i) != -INF && abs(master_lp_.lb->coeff(i) - res.x[i]) < 1e-8)
//      xlb_activity.coeffRef(i) = true;
//    if (master_lp_.ub && master_lp_.ub->coeff(i) != INF && abs(master_lp_.ub->coeff(i) - res.x[i]) < 1e-8)
//      xub_activity.coeffRef(i) = true;
//  }
//  Eigen::Matrix<bool, Eigen::Dynamic, 1> xlb_activity = res.dual_x_lb.array() > 0;
//  Eigen::Matrix<bool, Eigen::Dynamic, 1> xub_activity = res.dual_x_ub.array() > 0;

  // Count all previous active constraints at each index. This will be used in the next step to construct the complete
  // sets of active and inactive constraints
  Eigen::VectorXi active_count_ub(n_con_ub);

  int count_active = 0;
  for (size_t i = 0; i < n_con_ub; ++i) {
    active_count_ub.coeffRef(i) = count_active;
    if (ub_activity.coeff(i)) ++count_active;
  }

  long n_eq = (master_lp_.A_eq) ? master_lp_.A_eq->rows() : 0;
  long n_ub = ub_activity.count();
  long n_xlb = xlb_activity.count();
  long n_xub = xub_activity.count();
  n_active_ = n_eq + n_ub + n_xlb + n_xub;

  // Construct set of active constraints
  // TODO: simplify the following or make global
  Eigen::Matrix<bool, Eigen::Dynamic, 1> bone_eq(n_eq);
  bone_eq.setConstant(true);
  Eigen::SparseMatrix<double> eye(n_dim_master_, n_dim_master_);
  eye.setIdentity();
  Eigen::SparseMatrix<double> negative_eye = eye;
  negative_eye *= -1;

  Eigen::SparseMatrix<double> A_active = util::SparseOps::vstack_rows<double>(
      {
//        std::make_pair(*master_lp_.A_eq, bone_eq), // TODO: Need to check for existence before assignment.
//       std::make_pair(*master_lp_.A_ub, ub_activity), // TODO: Need to check for existence before assignment.
       std::make_pair(negative_eye, xlb_activity),
       std::make_pair(eye, xub_activity)}
  );

  // Add active constraints to residual matrix
  std::shared_ptr<Eigen::SparseMatrix<double>> A_activeT = std::make_shared<Eigen::SparseMatrix<double>>(
      A_active.transpose());
  A_activeT->conservativeResize(n_dim_master_ + 1, A_activeT->cols());
  std::shared_ptr<Eigen::VectorXd> lb = std::make_shared<Eigen::VectorXd>(n_active_);
  lb->setZero();
  for (size_t i = 0; i < n_eq; ++i) lb->coeffRef(i) = -INF;
  std::shared_ptr<Eigen::VectorXd> c = std::make_shared<Eigen::VectorXd>(n_active_);
  if (n_active_ > 0) residual_solver_->getQP().addVars(nullptr, nullptr, A_activeT, lb, nullptr);

  // Add variable for first visit
  std::shared_ptr<Eigen::VectorXd> a_eq = std::make_shared<Eigen::VectorXd>(n_dim_master_ + 1);
  a_eq->head(n_dim_master_) = cost_alpha_;
  a_eq->coeffRef(n_dim_master_) = 1.;

  residual_solver_->getQP().addVar(0, nullptr, a_eq, 0, INF);

  // Determine gradient adjustment
  Eigen::SparseMatrix<double> grad_adjust = Eigen::SparseMatrix<double>(n_dim_master_, n_dim_master_);
  grad_adjust.setIdentity();
  if (n_active_ > 0) {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> AAt;
    AAt.compute(A_active * A_active.transpose());
    grad_adjust -= A_active.transpose() * AAt.solve(A_active);
  }

  return grad_adjust * cost_alpha_;
}

const Eigen::VectorXd stuka::dLP::CRE::computeResidualSubsequent() {
  n_visit_++;

  std::shared_ptr<Eigen::VectorXd> a_eq = std::make_shared<Eigen::VectorXd>(n_dim_master_ + 1);
  a_eq->head(n_dim_master_) = cost_alpha_;
  a_eq->coeffRef(n_dim_master_) = 1.;

  residual_solver_->getQP().addVar(0., nullptr, a_eq, 0., INF);

  OptimizeState res;
  try {
    res = residual_solver_->solve();
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
