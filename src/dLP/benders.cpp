//
// Created by Avinash Madavan on 1/15/19.
//

#include <stuka/dLP/benders.h>

stuka::dLP::BendersDecomposition::BendersDecomposition(const stuka::dLP::DecomposedLinearProgram &dlp,
                                                       const stuka::Options &opts) : BaseDLPSolver(dlp, opts) {
  n_sub_calls_ = 0;
  n_sub_ = dlp.c.size() - 1;
  n_dim_master_ = dlp.c.back()->size();
  size_t n_con_ub_master = (dlp.b_ub.back()) ? dlp.b_ub.back()->size() : 0;
  size_t n_con_eq_master = (dlp.b_eq.back()) ? dlp.b_eq.back()->size() : 0;

  // Generate master problem
  master_lp_.c = std::make_shared<Eigen::VectorXd>(*dlp.c.back());
  master_lp_.c->conservativeResize(n_dim_master_ + n_sub_);
  for (size_t i = n_dim_master_; i < n_dim_master_ + n_sub_; ++i)
    master_lp_.c->coeffRef(i) = 1.;

  if (n_con_ub_master > 0) {
    master_lp_.A_ub = std::make_shared<Eigen::SparseMatrix<double>>(*dlp.A_ub.back());
    master_lp_.A_ub->conservativeResize(n_con_ub_master, n_dim_master_ + n_sub_);
    master_lp_.b_ub = dlp.b_ub.back();
  }

  if (n_con_eq_master > 0) {
    master_lp_.A_eq = std::make_shared<Eigen::SparseMatrix<double>>(*dlp.A_eq.back());
    master_lp_.A_eq->conservativeResize(n_con_eq_master, n_dim_master_ + n_sub_);
    master_lp_.b_eq = dlp.b_eq.back();
  }

  if (dlp.lb.back()) {
    master_lp_.lb = std::make_shared<Eigen::VectorXd>(*dlp.lb.back());
    master_lp_.lb->conservativeResize(n_dim_master_ + n_sub_);
    for (size_t i = n_dim_master_; i < n_dim_master_ + n_sub_; ++i)
      master_lp_.lb->coeffRef(i) = 0;
  }

  if (dlp.ub.back()) {
    master_lp_.ub = std::make_shared<Eigen::VectorXd>(*dlp.ub.back());
    master_lp_.ub->conservativeResize(n_dim_master_ + n_sub_);
    for (size_t i = n_dim_master_; i < n_dim_master_ + n_sub_; ++i)
      master_lp_.ub->coeffRef(i) = INF;
  }

  master_solver_ = util::createSolver(master_lp_, opts);

  // Set subproblem values
  subproblem_values_ = Eigen::VectorXd(n_sub_);
  subproblem_values_.setConstant(INF);

  // Set subproblems
  subproblems_.reserve(n_sub_);
  for (size_t i = 0; i < n_sub_; ++i) {
    Subproblem sub
        {dlp.c[i], dlp.A_ub[i], dlp.b_ub[i], dlp.C_ub[i], dlp.A_eq[i], dlp.b_eq[i], dlp.C_eq[i], dlp.lb[i], dlp.ub[i]};

    subproblems_.emplace_back(BendersSubproblem(std::move(sub), opts));
  }

  // Set initial point
  x_ = Eigen::VectorXd(n_dim_master_ + n_sub_);
  x_.setConstant(1.);
  if (opts.x0.size() == n_dim_master_) {
    x_.setConstant(0.);
    x_.head(n_dim_master_) = opts.x0;
  } else if (dlp.c.back() && dlp.c.back()->maxCoeff() > 1e-8 && dlp.c.back()->minCoeff() < -1e-8) {
    LP::LinearProgram tmp;
    tmp.c = dlp.c.back();
    tmp.A_ub = dlp.A_ub.back();
    tmp.b_ub = dlp.b_ub.back();
    tmp.A_eq = dlp.A_eq.back();
    tmp.b_eq = dlp.b_eq.back();
    tmp.lb = dlp.lb.back();
    tmp.ub = dlp.ub.back();

    OptimizeState res = util::createSolver(tmp, opts)->solve();
    if (res.status == 2) {
      x_.setConstant(0.);
      x_.head(n_dim_master_) = res.x;
    }
  }

}

void stuka::dLP::BendersDecomposition::iterate() {

  std::vector<BendersCut> cuts(n_sub_);
  Eigen::VectorXd x = x_.head(n_dim_master_);

  for (size_t i = 0; i < n_sub_; ++i) {
    n_sub_calls_++;
    cuts[i] = subproblems_[i].getBendersCut(x);
    subproblem_values_.coeffRef(i) = subproblems_[i].getValue();
  }

  OptimizeState res = solveMasterProblem(cuts);

  x_ = res.x;
}

bool stuka::dLP::BendersDecomposition::terminate() {
  // TODO: set tolerate from options
  return (x_.tail(n_sub_) - subproblem_values_).norm() < 1e-8;
}

const stuka::OptimizeState stuka::dLP::BendersDecomposition::getState() {
  stuka::OptimizeState res;

  res.x = x_.head(n_dim_master_);
  res.fun = master_lp_.c->dot(x_);
  res.error = (x_.tail(n_sub_) - subproblem_values_).norm();
  res.status = (res.error < 1e-8) ? 2 : 1;
  res.nit_sub = n_sub_calls_;

  return res;
}

stuka::OptimizeState
stuka::dLP::BendersDecomposition::solveMasterProblem(const std::vector<stuka::dLP::BendersCut> &cuts) {

  // Concatenate cuts and upper-bounding constraints
  std::shared_ptr<Eigen::SparseMatrix<double>> A_ub =
      std::make_shared<Eigen::SparseMatrix<double>>(n_sub_, n_dim_master_ + n_sub_);

  for (size_t i = 0; i < n_dim_master_; ++i) {
    A_ub->startVec(i);
    for (size_t j = 0; j < n_sub_; ++j)
      if (abs(cuts[j].a[i]) > 1e-16)
        A_ub->insertBack(j, i) = cuts[j].a.coeff(i);
  }

  for (size_t i = 0; i < n_sub_; ++i) {
    A_ub->startVec(n_dim_master_ + i);
    A_ub->insertBack(i, n_dim_master_ + i) = -1.;
  }

  A_ub->finalize();

  std::shared_ptr<Eigen::VectorXd> b_ub = std::make_shared<Eigen::VectorXd>(n_sub_);
  for (size_t i = 0; i < n_sub_; ++i)
    b_ub->coeffRef(i) = (abs(cuts[i].b) > 1e-8) ? cuts[i].b : 0;

  master_solver_->getLP().addConstrs_ub(A_ub, b_ub);

  OptimizeState res;
  try {
    res = master_solver_->solve();
  } catch (std::exception &e) {
    throw std::runtime_error("iterate: unable to solve master problem");
  }
  if (res.status != 2) {
    std::cout << "Failed to solve master problem: " << res.status << std::endl;
    throw std::runtime_error("iterate: unable to solve master problem");
  }

  return res;
}

