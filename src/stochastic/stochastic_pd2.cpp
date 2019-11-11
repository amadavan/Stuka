//
// Created by Avinash Madavan on 2019-06-26.
//

//
// Created by Avinash Madavan on 2019-05-07.
//

#include <stuka/stochastic/stochastic_pd2.h>

stuka::stochastic::StochasticPrimalDual2::StochasticPrimalDual2(const stuka::stochastic::Program &prog,
                                                                const stuka::Options &opts)
    : BaseStochasticSolver(prog, opts), prog_(prog), measure_(opts.measure) {
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

  x_ = measure_->augmentState(x_);

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

  double step = 10. / sqrt(nit_ + 0.);

  Eigen::VectorXd wk = prog_.sample();
  x_ -= step_ * measure_->df(prog_, x_, wk);
  if (prog_.g) x_ -= step_ * measure_->dg(prog_, x_, wk) * z_.head(m_g_);
  if (prog_.h) x_ -= step_ * prog_.dh(x_) * z_.tail(m_h_);
  x_ = measure_->projX(prog_, x_);

  wk = prog_.sample();
  if (prog_.g) z_.head(m_g_) += step_ * measure_->g(prog_, x_, wk);
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
