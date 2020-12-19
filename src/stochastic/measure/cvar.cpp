//
// Created by Avinash Madavan on 10/31/19.
//

#include <stuka/stochastic/measure/cvar.h>

double stuka::stochastic::ConditionalValueAtRisk::f(const Program &p,
                                                    const Eigen::VectorXd &xu,
                                                    const Eigen::VectorXd &xi) {
  Eigen::VectorXd x = xu.head(n_);
  double u0 = xu.coeff(n_);
  Eigen::VectorXd u = xu.tail(m_g_);

  return u0 + alphap_ * std::max(0., p.f(x, xi) - u0);
}

Eigen::VectorXd stuka::stochastic::ConditionalValueAtRisk::g(const Program &p,
                                                             const Eigen::VectorXd &xu,
                                                             const Eigen::VectorXd &xi) {
  Eigen::VectorXd x = xu.head(n_);
  double u0 = xu.coeff(n_);
  Eigen::VectorXd u = xu.tail(m_g_);

  Eigen::VectorXd gmu = p.g(x, xi) - u;
  gmu = (gmu.array() < 0).select(0, gmu);

  return u + betap_.cwiseProduct(gmu);
}

Eigen::VectorXd stuka::stochastic::ConditionalValueAtRisk::df(const Program &p,
                                                              const Eigen::VectorXd &xu,
                                                              const Eigen::VectorXd &xi) {
  Eigen::VectorXd x = xu.head(n_);
  double u0 = xu.coeff(n_);
  Eigen::VectorXd u = xu.tail(m_g_);

  Eigen::VectorXd df = Eigen::VectorXd::Zero(xu.size());
  df.coeffRef(n_) = 1.;

  if (p.f(x, xi) >= u0) {
    df.head(n_) = alphap_ * p.df(x, xi);
    df.coeffRef(n_) -= alphap_;
  }

  return df;
}

Eigen::MatrixXd stuka::stochastic::ConditionalValueAtRisk::dg(const Program &p,
                                                              const Eigen::VectorXd &xu,
                                                              const Eigen::VectorXd &xi) {
  Eigen::VectorXd x = xu.head(n_);
  double u0 = xu.coeff(n_);
  Eigen::VectorXd u = xu.tail(m_g_);

  Eigen::MatrixXd dg = Eigen::MatrixXd::Zero(xu.size(), m_g_);
  Eigen::MatrixXd prog_dg = p.dg(x, xi);

  Eigen::Array<bool, Eigen::Dynamic, 1> gpos((p.g(x, xi) - u).array() >= 0.);
  for (size_t i = 0; i < m_g_; ++i) {
    dg.coeffRef(n_ + 1 + i, i) = 1.;
    if (gpos.coeff(i)) {
      dg.block(0, i, n_, 1) += betap_.coeff(i) * prog_dg.col(i);
      dg.coeffRef(n_ + 1 + i, i) -= betap_.coeff(i);
    }
  }

  return dg;
}

Eigen::VectorXd stuka::stochastic::ConditionalValueAtRisk::projX(const Program &p, const Eigen::VectorXd &xu) {
  Eigen::VectorXd xu_proj(xu);
  Eigen::VectorXd x = xu.head(n_);
  double u0 = xu.coeff(n_);
  Eigen::VectorXd u = xu.tail(m_g_);

  xu_proj.head(n_) = p.projX(x);
  for (unsigned long i = 0; i < u.size(); ++i) {
    if (u[i] < -g_max_.coeff(i)) xu_proj.coeffRef(n_ + 1 + i) = -g_max_.coeff(i);
    if (u[i] > g_max_.coeff(i)) xu_proj.coeffRef(n_ + 1 + i) = g_max_.coeff(i);
  }

  return xu_proj;
}

Eigen::VectorXd stuka::stochastic::ConditionalValueAtRisk::augmentState(Eigen::VectorXd &x) {
  x.conservativeResize(n_ + 1 + m_g_);
  for (size_t i = 0; i < 1 + m_g_; ++i) x.coeffRef(i) = 0.;
  return x;
}
