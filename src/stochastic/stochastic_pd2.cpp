//
// Created by Avinash Madavan on 2019-06-26.
//

//
// Created by Avinash Madavan on 2019-05-07.
//

#include <stuka/stochastic/stochastic_pd2.h>

stuka::stochastic::StochasticPrimalDual2::StochasticPrimalDual2(const stuka::stochastic::Program &prog,
                                                              const stuka::Options &opts) : BaseStochasticSolver(prog,
                                                                                                                 opts),
                                                                                            prog_(prog) {
  step_ = (opts.step != 0) ? opts.step : 1.;

  if (opts.x0.size() > 0)
    x_ = opts.x0;
  else
    throw std::runtime_error("StochasticPrimalDual2::StochasticPrimalDual: Initial point required.");

  // Compute initial sample for sizing
  Eigen::VectorXd xi = prog_.sample();

  n_ = x_.size();
  m_g_ = (prog_.g) ? prog_.g(x_, xi).rows() : 0;
  m_h_ = (prog_.h) ? prog_.h(x_).rows() : 0;

  z_ = Eigen::VectorXd(m_g_ + m_h_);
  z_.setZero();

  if (prog_.alpha == 0. && (prog_.beta.size() == 0. || prog_.beta.maxCoeff() < 1e-16)) {
    f_ = prog_.f;
    g_ = prog_.g;
    df_ = prog_.df;
    dg_ = prog_.dg;
    projX_ = prog_.projX;
  } else {
    x_.conservativeResize(n_ + 1 + m_g_);
    for (size_t i = 0; i < 1 + m_g_; ++i) x_.coeffRef(i) = 0.;
    alphap_ = 1. / (1. - prog_.alpha);
    betap_ = 1. / (1. - prog_.beta.array());
    f_ = std::bind(&StochasticPrimalDual2::f_cvar, this, std::placeholders::_1, std::placeholders::_2);
    g_ = std::bind(&StochasticPrimalDual2::g_cvar, this, std::placeholders::_1, std::placeholders::_2);
    df_ = std::bind(&StochasticPrimalDual2::df_cvar, this, std::placeholders::_1, std::placeholders::_2);
    dg_ = std::bind(&StochasticPrimalDual2::dg_cvar, this, std::placeholders::_1, std::placeholders::_2);
    projX_ = std::bind(&StochasticPrimalDual2::projX_cvar, this, std::placeholders::_1);
  }

  xbar_ = x_;
  zbar_ = z_;
  nit_ = 1;

  xbarp_ = x_;
  xbarpp_ = xbarp_;
  xbarppp_ = xbarpp_;

  stepsum_ = 0.;

  restart_ = 0;
}

void stuka::stochastic::StochasticPrimalDual2::iterate() {

  xbarppp_ = xbarpp_;
  xbarpp_ = xbarp_;
  xbarp_ = xbar_;

  Eigen::VectorXd wk = prog_.sample();
  double step = 10./sqrt(nit_ + 0.);

  x_ -= step * df_(x_, wk);
  if (prog_.g) x_ -= step * dg_(x_, wk) * z_.head(m_g_);
  if (prog_.h) x_ -= step * prog_.dh(x_) * z_.tail(m_h_);
  x_ = projX_(x_);

  wk = prog_.sample();
  if (prog_.g) z_.head(m_g_) += step * g_(x_, wk);
  if (prog_.h) z_.tail(m_h_) += step * prog_.h(x_);
  z_ = (z_.array() < 0).select(0, z_);

  xbar_ = (xbar_ * stepsum_ + x_ * step) / (stepsum_ + step);
  zbar_ = (zbar_ * stepsum_ + z_ * step) / (stepsum_ + step);
  stepsum_ += step;
  nit_++;

}

bool stuka::stochastic::StochasticPrimalDual2::terminate() {
  if ((xbarp_ - xbar_).norm() < 1e-4)
    if ((xbarpp_ - xbar_).norm() < 1e-4)
      if ((xbarppp_ - xbar_).norm() < 1e-4) {
        if (restart_ == 3)
          return true;
        else {
          nit_ = 1;
          restart_++;
          return false;
        }
      }
  return false;
}

const stuka::OptimizeState stuka::stochastic::StochasticPrimalDual2::getState() {
  return {
      xbar_,
      zbar_
  };
}

double stuka::stochastic::StochasticPrimalDual2::f_cvar(Eigen::VectorXd &xu, Eigen::VectorXd &xi) {
  Eigen::VectorXd x = xu.head(n_);
  double u0 = xu.coeff(n_);
  Eigen::VectorXd u = xu.tail(m_g_);

  return u0 + (1. / (1. - alphap_)) * std::max(0., prog_.f(x, xi) - u0);
}

Eigen::VectorXd stuka::stochastic::StochasticPrimalDual2::g_cvar(Eigen::VectorXd &xu, Eigen::VectorXd &xi) {
  Eigen::VectorXd x = xu.head(n_);
  double u0 = xu.coeff(n_);
  Eigen::VectorXd u = xu.tail(m_g_);

  Eigen::VectorXd gmu = prog_.g(x, xi) - u;
  gmu = (gmu.array() < 0).select(0, gmu);

  return u + betap_.cwiseProduct(gmu);
}

Eigen::VectorXd stuka::stochastic::StochasticPrimalDual2::df_cvar(Eigen::VectorXd &xu, Eigen::VectorXd &xi) {
  Eigen::VectorXd x = xu.head(n_);
  double u0 = xu.coeff(n_);
  Eigen::VectorXd u = xu.tail(m_g_);

  Eigen::VectorXd df = Eigen::VectorXd(xu.size());
  df.setZero();
  df.coeffRef(n_) = 1.;

  if (prog_.f(x, xi) >= u0) {
    df.head(n_) = alphap_ * prog_.df(x, xi);
    df.coeffRef(n_) -= alphap_;
  }

  return df;
}

Eigen::MatrixXd stuka::stochastic::StochasticPrimalDual2::dg_cvar(Eigen::VectorXd &xu, Eigen::VectorXd &xi) {
  Eigen::VectorXd x = xu.head(n_);
  double u0 = xu.coeff(n_);
  Eigen::VectorXd u = xu.tail(m_g_);

  Eigen::MatrixXd dg = Eigen::MatrixXd(xu.size(), m_g_);
  Eigen::MatrixXd prog_dg = prog_.dg(x, xi);
  dg.setZero();

  Eigen::Matrix<bool, Eigen::Dynamic, 1> gpos((prog_.g(x, xi) - u).array() >= 0.);
  for (size_t i = 0; i < m_g_; ++i) {
    dg.coeffRef(n_ + 1 + i, i) = 1.;
    if (gpos.coeff(i)) {
      dg.block(0, i, n_, 1) += betap_.coeff(i) * prog_dg.col(i);
      dg.coeffRef(n_ + 1 + i, i) -= betap_.coeff(i);
    }
  }

  return dg;
}

Eigen::VectorXd stuka::stochastic::StochasticPrimalDual2::projX_cvar(Eigen::VectorXd &xu) {
  Eigen::VectorXd x = xu.head(n_);
  double u0 = xu.coeff(n_);
  Eigen::VectorXd u = xu.tail(m_g_);

  xu.head(n_) = prog_.projX(x);
  for (unsigned long i = 0; i < u.size(); ++i) {
    if (u[i] < -prog_.gmax.coeff(i)) xu.coeffRef(n_ + 1 + i) = -prog_.gmax.coeff(i);
    if (u[i] > prog_.gmax.coeff(i)) xu.coeffRef(n_ + 1 + i) = prog_.gmax.coeff(i);
  }

  return xu;
}
