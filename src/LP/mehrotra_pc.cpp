//
// Created by Avinash Madavan on 2019-03-02.
//

#include <stuka/LP/mehrotra_pc.h>

stuka::LP::MehrotraPC::MehrotraPC(const stuka::LP::LinearProgram &lp, const stuka::Options &opts)
    : BaseLPSolver(lp, opts), prog_(lp), eps_(opts.tol) {

  // Useful vectors
  double bnorm = prog_.b_->norm(), cnorm = prog_.c_->norm();
  bc_ = (bnorm > cnorm) ? bnorm + 1. : cnorm + 1.;

  e_ = Eigen::VectorXd(prog_.n_dim_);
  e_.setOnes();

  // Construct AAt
  Eigen::SparseMatrix<double> M = *prog_.A_ * prog_.A_->transpose();
  solver_.analyzePattern(M);
  solver_.factorize(M);

  // min norm(x) s.t. Ax = b
  x_ = prog_.A_->transpose() * solver_.solve(*prog_.b_);

  // min norm(s) s.t. A'y + s = c
  y_ = solver_.solve((*prog_.A_) * (*prog_.c_));
  s_ = *prog_.c_ - prog_.A_->transpose() * y_;

  double dx = -1.5 * x_.minCoeff();
  double ds = -1.5 * s_.minCoeff();

  // Ensure positivity
  dx = (dx > 0) ? dx : 0;
  ds = (ds > 0) ? ds : 0;

  double pdct = 0.5 * (x_ + dx * e_).dot(s_ + ds * e_);
  double dx_c = dx + pdct / (s_.sum() + prog_.n_dim_ * ds);
  double ds_c = ds + pdct / (x_.sum() + prog_.n_dim_ * dx);

  // Compute initial point
  x_ += dx_c * e_;
  s_ += ds_c * e_;
}

void stuka::LP::MehrotraPC::iterate() {

  Eigen::VectorXd dx, dy, ds;
  double alpha_p, alpha_d;

  mu_ = x_.dot(s_) / prog_.n_dim_;
  Rc_ = *prog_.c_ - prog_.A_->transpose() * y_ - s_;
  Rb_ = *prog_.b_ - *prog_.A_ * x_;
  Rxs_ = -x_.cwiseProduct(s_);

  // Compute the affine predictor step
  std::tie(dx, dy, ds) = solveKKT(Rc_, Rb_, Rxs_, true);
  std::tie(alpha_p, alpha_d) = computeStepSize(dx, ds, 1.);

  double mu_aff = 1. / prog_.n_dim_ * (x_ + alpha_p * dx).dot(s_ + alpha_d * ds);
  double sigma = pow(mu_aff / mu_, 3);

  // Compute the corrector step
  std::tie(dx, dy, ds) = solveKKT(Rc_, Rb_, Rxs_ + sigma * mu_ * e_ - dx.cwiseProduct(ds), false);
  std::tie(alpha_p, alpha_d) = computeStepSize(dx, ds, 0.99995);

  x_ += alpha_p * dx;
  y_ += alpha_d * dy;
  s_ += alpha_d * ds;

}

bool stuka::LP::MehrotraPC::terminate() {
  error_ = (Rc_.norm() + Rb_.norm() + Rxs_.norm()) / bc_;
  return error_ < eps_;
}

const stuka::OptimizeState stuka::LP::MehrotraPC::getState() {
  OptimizeState res;
  res.x = prog_.revertState(x_);
  res.dual_eq = y_;
  res.dual_ub = s_;
  res.error = error_;

  return res;
}

stuka::LP::BaseLinearProgram &stuka::LP::MehrotraPC::getLP() {
  throw std::runtime_error("LP::MehrotraPC currently does not support modification of LP");
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>
stuka::LP::MehrotraPC::solveKKT(Eigen::VectorXd Rc, Eigen::VectorXd Rb, Eigen::VectorXd Rxs, bool resolve) {
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> sdiag_inv = s_.asDiagonal().inverse();
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> xdiag = x_.asDiagonal();

  if (resolve) {
    solver_.factorize((*prog_.A_) * sdiag_inv * xdiag * prog_.A_->transpose());
  }

  Eigen::VectorXd rhs_y = *prog_.A_ * sdiag_inv * xdiag * Rc + Rb - *prog_.A_ * sdiag_inv * Rxs;

  Eigen::VectorXd dy = solver_.solve(rhs_y);
  Eigen::VectorXd ds = Rc - prog_.A_->transpose() * dy;
  Eigen::VectorXd dx = sdiag_inv * (Rxs - x_.cwiseProduct(ds));

  return std::make_tuple(dx, dy, ds);
}

std::tuple<double, double>
stuka::LP::MehrotraPC::computeStepSize(Eigen::VectorXd dx, Eigen::VectorXd ds, double scale) {
  // Primal step
  double alpha_p = dx.cwiseQuotient(x_).minCoeff();
  alpha_p = (alpha_p > -1.) ? scale : -scale / alpha_p;
  alpha_p = (alpha_p > 1.) ? 1. : alpha_p;

  // Dual step
  double alpha_d = ds.cwiseQuotient(s_).minCoeff();
  alpha_d = (alpha_d > -1.) ? scale : -scale / alpha_d;
  alpha_d = (alpha_d > 1.) ? 1. : alpha_d;

  return std::make_tuple(alpha_p, alpha_d);
}
