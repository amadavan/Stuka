//
// Created by Avinash Madavan on 10/8/19.
//

#ifndef STUKA_STOCHASTIC_MEASURE_G_ENTROPIC_H_
#define STUKA_STOCHASTIC_MEASURE_G_ENTROPIC_H_

// TODO: find a better way to handle different divergences/functions for each constraint

#include <functional>

#include "measure.h"

namespace stuka { namespace stochastic {

class gConjugateFunction {
 public:
  virtual double g(const double &x) = 0;
  virtual double dg(const double &x) = 0;
};

template<class gFunction>
class gEntropicRisk : public MeasurableProgram {
 public:
  gEntropicRisk(const Program &prog, const Eigen::VectorXd &divergence)
      : MeasurableProgram(prog) {

    gFunction func = gFunction();

    f = [this, prog, func, divergence](const Eigen::VectorXd &x_ent, const Eigen::VectorXd &xi) {
      Eigen::VectorXd x = x_ent.head(n_);
      double u0 = x_ent.coeff(n_);
      double t0 = x_ent.coeff(n_ + m_g_);

      return u0 + t0 * func.g((prog.f(x, xi) - u0) / t0 + divergence);
    };

    g = [this, prog, func, divergence](const Eigen::VectorXd &x_ent, const Eigen::VectorXd &xi) {
      Eigen::VectorXd x = x_ent.head(n_);
      Eigen::VectorXd u = x_ent.segment(n_, m_g_);
      Eigen::VectorXd t = x_ent.tail(m_g_);

      Eigen::VectorXd g_val = u;
      Eigen::VectorXd g_ent = (prog.g(x, xi) - u) / t + divergence;

      for (size_t i = 0; i < m_g_; ++i) {
        g_val.coeffRef(i) += t * func.g(g_ent.coeff(i))
      }

      return g_val;
    };

    df = [this, prog, func, divergence](const Eigen::VectorXd &x_ent, const Eigen::VectorXd &xi) {
      Eigen::VectorXd x = x_ent.head(n_);
      double u0 = x_ent.coeff(n_);
      double t0 = x_ent.coeff(n_ + 1 + m_g_);

      Eigen::VectorXd df = Eigen::VectorXd(x_ent.size());
      df.setZero();

      double x_g = prog.f(x, xi) - u0)/t0;
      double g_eval = func.dg(x_g + divergence.coeff(0));

      df.head(n_) = g_eval * prog.df(x, xi);
      df.coeffRef(n_) = 1 - g_eval;
      df.coeffRef(n_ + m_g_) = g_eval - g_eval * x_g;

      return df;
    };

    dg = [this, prog, func, divergence](const Eigen::VectorXd &x_ent, const Eigen::VectorXd &xi) {
      Eigen::VectorXd x = x_ent.head(n_);
      Eigen::VectorXd u = x_ent.segment(n_ + 1, m_g_);
      Eigen::VectorXd t = x_ent.tail(m_g_);

      Eigen::MatrixXd dg = Eigen::MatrixXd(x_ent.size(), m_g_);
      Eigen::MatrixXd prog_dg = prog.dg(x, xi);
      dg.setZero();

      Eigen::VectorXd x_g = prog.g(x, xi) - u)/t;

      for (size_t i = 0; i < m_g_; ++i) {
        double g_eval = func.g(x_g.coeff(i) + divergence.coeff(i + 1));

        dg.block(0, i, n_, 1) = g_eval * prog_dg.col(i);
        dg.coeffRef(n_ + 1 + i, i) = 1 - g_eval;
        dg.coeffRef(n_ + 1 + m_g_ + i, i) = g_eval - g_eval * x_g.coeff(i);
      }

      return dg;
    };

    projX = [this, prog](const Eigen::VectorXd &x_ent) {
      Eigen::VectorXd x = x_ent.head(n_);

      Eigen::VectorXd x_proj(x_ent);
      x_proj.head(n_) = prog.projX(x_ent);

      return x_proj;
    }
  };

 private:
  size_t n_;
  size_t m_g_;
};

}}

#endif //STUKA_STOCHASTIC_MEASURE_G_ENTROPIC_H_
