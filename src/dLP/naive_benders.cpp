//
// Created by Avinash on 8/31/2020.
//

//
// Created by Avinash Madavan on 1/15/19.
//

#include <stuka/dLP/naive_benders.h>

stuka::dLP::NaiveBendersDecomposition::NaiveBendersDecomposition(const stuka::dLP::DecomposedLinearProgram &dlp,
                                                       const stuka::Options &opts) : BaseDLPSolver(dlp, opts),
                                                                                     opts_(opts) {
  n_sub_calls_ = 0;
  n_sub_ = dlp.c.size() - 1;
  n_dim_master_ = dlp.c.back()->size();
  size_t n_con_ub_master = (dlp.b_ub.back()) ? dlp.b_ub.back()->size() : 0;
  size_t n_con_eq_master = (dlp.b_eq.back()) ? dlp.b_eq.back()->size() : 0;

  // Generate master problem
  master_lp_.c = util::DenseOps::unique_copy(dlp.c.back());
  master_lp_.c->conservativeResize(n_dim_master_ + n_sub_);
  for (size_t i = n_dim_master_; i < n_dim_master_ + n_sub_; ++i)
    master_lp_.c->coeffRef(i) = 1.;

  if (n_con_ub_master > 0) {
    master_lp_.A_ub = util::SparseOps::unique_copy(dlp.A_ub.back());
    master_lp_.A_ub->conservativeResize(n_con_ub_master, n_dim_master_ + n_sub_);
    master_lp_.b_ub = util::DenseOps::unique_copy(dlp.b_ub.back());
  }

  if (n_con_eq_master > 0) {
    master_lp_.A_eq = util::SparseOps::unique_copy(dlp.A_eq.back());
    master_lp_.A_eq->conservativeResize(n_con_eq_master, n_dim_master_ + n_sub_);
    master_lp_.b_eq = util::DenseOps::unique_copy(dlp.b_eq.back());
  }

  if (dlp.lb.back()) {
    master_lp_.lb = util::DenseOps::unique_copy(dlp.lb.back());
    master_lp_.lb->conservativeResize(n_dim_master_ + n_sub_);
    for (size_t i = n_dim_master_; i < n_dim_master_ + n_sub_; ++i)
      master_lp_.lb->coeffRef(i) = 0;
  }

  if (dlp.ub.back()) {
    master_lp_.ub = util::DenseOps::unique_copy(dlp.ub.back());
    master_lp_.ub->conservativeResize(n_dim_master_ + n_sub_);
    for (size_t i = n_dim_master_; i < n_dim_master_ + n_sub_; ++i)
      master_lp_.ub->coeffRef(i) = INF;
  }

  // Set subproblem values
  subproblem_values_ = Eigen::VectorXd::Constant(n_sub_, INF);

  // Set subproblems
  subproblems_.reserve(n_sub_);
  Options opts_sub(opts_);
  opts_sub.lazy = opts_.lazy_sub;

  for (size_t i = 0; i < n_sub_; ++i) {
    subproblems_.emplace_back((Subproblem){
        util::DenseOps::unique_copy(dlp.c[i]),
        util::SparseOps::unique_copy(dlp.A_ub[i]),
        util::DenseOps::unique_copy(dlp.b_ub[i]),
        util::SparseOps::unique_copy(dlp.C_ub[i]),
        util::SparseOps::unique_copy(dlp.A_eq[i]),
        util::DenseOps::unique_copy(dlp.b_eq[i]),
        util::SparseOps::unique_copy(dlp.C_eq[i]),
        util::DenseOps::unique_copy(dlp.lb[i]),
        util::DenseOps::unique_copy(dlp.ub[i])
    }, opts_sub);
  }

  // Set initial point
  x_ = Eigen::VectorXd::Constant(n_dim_master_ + n_sub_, 1);

  if (opts_.x0.size() == n_dim_master_) {
    x_.setConstant(0.);
    x_.head(n_dim_master_) = opts_.x0;
  } else if (dlp.c.back() && dlp.c.back()->maxCoeff() > 1e-8 && dlp.c.back()->minCoeff() < -1e-8) {
    LP::LinearProgram tmp;
    tmp.c = util::DenseOps::unique_copy(dlp.c.back());
    tmp.A_ub = util::SparseOps::unique_copy(dlp.A_ub.back());
    tmp.b_ub = util::DenseOps::unique_copy(dlp.b_ub.back());
    tmp.A_eq = util::SparseOps::unique_copy(dlp.A_eq.back());
    tmp.b_eq = util::DenseOps::unique_copy(dlp.b_eq.back());
    tmp.lb = util::DenseOps::unique_copy(dlp.lb.back());
    tmp.ub = util::DenseOps::unique_copy(dlp.ub.back());

    OptimizeState res = util::createSolver(tmp, opts_)->solve();
    if (res.status == 2) {
      x_.setConstant(0.);
      x_.head(n_dim_master_) = res.x;
    }
  }

}

void stuka::dLP::NaiveBendersDecomposition::iterate() {

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

bool stuka::dLP::NaiveBendersDecomposition::terminate() {
  return (x_.tail(n_sub_) - subproblem_values_).norm() < opts_.tol;
}

const stuka::OptimizeState stuka::dLP::NaiveBendersDecomposition::getState() {
  stuka::OptimizeState res;

  res.x = x_.head(n_dim_master_);
  res.fun = master_lp_.c->dot(x_);
  res.error = (x_.tail(n_sub_) - subproblem_values_).norm();
  res.status = (res.error < 1e-8) ? 2 : 1;
  res.nit_sub = n_sub_calls_;

  return res;
}

stuka::OptimizeState
stuka::dLP::NaiveBendersDecomposition::solveMasterProblem(const std::vector<stuka::dLP::BendersCut> &cuts) {

  // Concatenate cuts and upper-bounding constraints
  std::unique_ptr<Eigen::SparseMatrix<double>> A_ub =
      std::make_unique<Eigen::SparseMatrix<double>>(n_sub_, n_dim_master_ + n_sub_);

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

  std::unique_ptr<Eigen::VectorXd> b_ub = std::make_unique<Eigen::VectorXd>(n_sub_);
  for (size_t i = 0; i < n_sub_; ++i)
    b_ub->coeffRef(i) = (abs(cuts[i].b) > 1e-8) ? cuts[i].b : 0;

  if (master_lp_.A_ub) {
    Eigen::SparseMatrix<double> _A_ub = util::SparseOps::vstack<double>({*master_lp_.A_ub, *A_ub});
    Eigen::VectorXd _b_ub = Eigen::VectorXd(master_lp_.b_ub->size() + b_ub->size());
    _b_ub.head(master_lp_.b_ub->size()) = *master_lp_.b_ub;
    _b_ub.tail(b_ub->size()) = *b_ub;

    master_lp_.A_ub = std::make_unique<Eigen::SparseMatrix<double>>(_A_ub);
    master_lp_.b_ub = std::make_unique<Eigen::VectorXd>(_b_ub);
  }
  else {
    master_lp_.A_ub = std::move(A_ub);
    master_lp_.b_ub = std::move(b_ub);
  }

  std::unique_ptr<LP::BaseLPSolver> master_solver = util::createSolver(master_lp_, opts_);

  OptimizeState res;
  try {
    res = master_solver->solve();
  } catch (std::exception &e) {
    throw std::runtime_error("iterate: unable to solve master problem");
  }
  if (res.status != 2) {
    std::cout << "Failed to solve master problem: " << res.status << std::endl;
    throw std::runtime_error("iterate: unable to solve master problem");
  }

  return res;
}

