//
// Created by Avinash Madavan on 2019-09-16.
//

#ifndef STUKA_STOCHASTIC_MEASURE_CVAR_H_
#define STUKA_STOCHASTIC_MEASURE_CVAR_H_

#include "measure.h"

namespace stuka { namespace stochastic {
class ConditionalValueAtRisk : public MeasurableProgram {
 public:
  ConditionalValueAtRisk(const Program &prog)
      : MeasurableProgram(prog), alphap_(1. / (1. - prog.alpha)), betap_(1. / (1. - prog.beta.array())) {
    m_g_ = betap_.size();

    f = [this, prog](const Eigen::VectorXd &xu, const Eigen::VectorXd &xi) {
      Eigen::VectorXd x = xu.head(n_);
      double u0 = xu.coeff(n_);
      Eigen::VectorXd u = xu.tail(m_g_);

      return u0 + alphap_ * std::max(0., prog.f(x, xi) - u0);
    };

    g = [this, prog](const Eigen::VectorXd &xu, const Eigen::VectorXd &xi) {
      Eigen::VectorXd x = xu.head(n_);
      double u0 = xu.coeff(n_);
      Eigen::VectorXd u = xu.tail(m_g_);

      Eigen::VectorXd gmu = prog.g(x, xi) - u;
      gmu = (gmu.array() < 0).select(0, gmu);

      return u + betap_.cwiseProduct(gmu);
    };

    df = [this, prog](const Eigen::VectorXd &xu, const Eigen::VectorXd &xi) {
      Eigen::VectorXd x = xu.head(n_);
      double u0 = xu.coeff(n_);
      Eigen::VectorXd u = xu.tail(m_g_);

      Eigen::VectorXd df = Eigen::VectorXd(xu.size());
      df.setZero();
      df.coeffRef(n_) = 1.;

      if (prog.f(x, xi) >= u0) {
        df.head(n_) = alphap_ * prog.df(x, xi);
        df.coeffRef(n_) -= alphap_;
      }

      return df;
    };

    dg = [this, prog](const Eigen::VectorXd &xu, const Eigen::VectorXd &xi) {
      Eigen::VectorXd x = xu.head(n_);
      double u0 = xu.coeff(n_);
      Eigen::VectorXd u = xu.tail(m_g_);

      Eigen::MatrixXd dg = Eigen::MatrixXd(xu.size(), m_g_);
      Eigen::MatrixXd progdg = prog.dg(x, xi);
      dg.setZero();

      Eigen::Matrix<bool, Eigen::Dynamic, 1> gpos((prog.g(x, xi) - u).array() >= 0.);
      for (size_t i = 0; i < m_g_; ++i) {
        dg.coeffRef(n_ + 1 + i, i) = 1.;
        if (gpos.coeff(i)) {
          dg.block(0, i, n_, 1) += betap_.coeff(i) * progdg.col(i);
          dg.coeffRef(n_ + 1 + i, i) -= betap_.coeff(i);
        }
      }

      return dg;
    };

    projX = [this, prog](const Eigen::VectorXd &xu) {
      Eigen::VectorXd xu_proj(xu);
      Eigen::VectorXd x = xu.head(n_);
      double u0 = xu.coeff(n_);
      Eigen::VectorXd u = xu.tail(m_g_);

      xu_proj.head(n_) = prog.projX(x);
      for (unsigned long i = 0; i < u.size(); ++i) {
        if (u[i] < -prog.gmax.coeff(i)) xu_proj.coeffRef(n_ + 1 + i) = -prog.gmax.coeff(i);
        if (u[i] > prog.gmax.coeff(i)) xu_proj.coeffRef(n_ + 1 + i) = prog.gmax.coeff(i);
      }

      return xu_proj;
    };
  }

 private:
  size_t n_;
  size_t m_g_;
  const size_t alphap_;
  const Eigen::VectorXd betap_;
};
}}

#endif //STUKA_STOCHASTIC_MEASURE_CVAR_H_
