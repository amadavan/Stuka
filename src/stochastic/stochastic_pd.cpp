//
// Created by Avinash Madavan on 2019-05-07.
//

#include <stuka/stochastic/stochastic_pd.h>

stuka::stochastic::StochasticPrimalDual::StochasticPrimalDual(const stuka::stochastic::Program &prog,
                                                              const stuka::Options &opts)
    : BaseStochasticSolver(prog, opts), prog_(prog), measure_(opts.measure) {
  step_ = (opts.step != 0) ? opts.step : 1.;

  if (opts.x0.size() > 0)
    x_ = opts.x0;
  else
    throw std::runtime_error("StochasticPrimalDual::StochasticPrimalDual: Initial point required.");

  // Compute initial sample for sizing
  Eigen::VectorXd xi = prog_.sample();

  n_ = x_.size();
  m_g_ = (prog_.g) ? prog_.g(x_, xi).rows() : 0;
  m_h_ = (prog_.h) ? prog_.h(x_).rows() : 0;

  z_ = Eigen::VectorXd::Zero(m_g_ + m_h_);

  x_ = measure_->augmentState(x_);

  x_bar_ = x_;
  z_bar_ = z_;
  nit_ = 1;

}

void stuka::stochastic::StochasticPrimalDual::iterate() {
  Eigen::VectorXd wk = prog_.sample();
  x_ -= step_ * measure_->df(prog_, x_, wk);
  if (prog_.g) x_ -= step_ * measure_->dg(prog_, x_, wk) * z_.head(m_g_);
  if (prog_.h) x_ -= step_ * prog_.dh(x_) * z_.tail(m_h_);
  x_ = measure_->projX(prog_, x_);

  wk = prog_.sample();
  if (prog_.g) z_.head(m_g_) += step_ * measure_->g(prog_, x_, wk);
  if (prog_.h) z_.tail(m_h_) += step_ * prog_.h(x_);
  z_ = (z_.array() < 0).select(0, z_);

  x_bar_ = (x_bar_ * nit_ + x_) / (nit_ + 1);
  z_bar_ = (z_bar_ * nit_ + z_) / (nit_ + 1);
  nit_++;

}

bool stuka::stochastic::StochasticPrimalDual::terminate() {
  return false;
}

const stuka::OptimizeState stuka::stochastic::StochasticPrimalDual::getState() {
  return {
      x_bar_,
      z_bar_
  };
}
