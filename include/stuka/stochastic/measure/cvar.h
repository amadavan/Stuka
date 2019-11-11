//
// Created by Avinash Madavan on 2019-09-16.
//

#ifndef STUKA_STOCHASTIC_MEASURE_CVAR_H_
#define STUKA_STOCHASTIC_MEASURE_CVAR_H_

#include "../program.h"
#include "measure.h"

namespace stuka { namespace stochastic {
class ConditionalValueAtRisk : public ProgramMeasure {
 public:
  ConditionalValueAtRisk(const size_t &n_var,
                         const double &alpha,
                         const Eigen::VectorXd &beta,
                         const Eigen::VectorXd &g_max)
      : n_(n_var),
        m_g_(g_max.size()),
        alphap_(1. / (1. - alpha)),
        betap_((beta.size() > 0) ? 1. / (1. - beta.array()) : Eigen::VectorXd(m_g_)),
        g_max_(g_max) {}
  ~ConditionalValueAtRisk() = default;

  double f(const Program &, const Eigen::VectorXd &, const Eigen::VectorXd &) override;
  Eigen::VectorXd g(const Program &, const Eigen::VectorXd &, const Eigen::VectorXd &) override;
  Eigen::VectorXd df(const Program &, const Eigen::VectorXd &, const Eigen::VectorXd &) override;
  Eigen::MatrixXd dg(const Program &, const Eigen::VectorXd &, const Eigen::VectorXd &) override;
  Eigen::VectorXd projX(const Program &, const Eigen::VectorXd &) override;
  Eigen::VectorXd augmentState(Eigen::VectorXd &x) override;

 private:
  const size_t n_;
  const size_t m_g_;
  const size_t alphap_;
  const Eigen::VectorXd betap_;
  const Eigen::VectorXd g_max_;
};
}}

#endif //STUKA_STOCHASTIC_MEASURE_CVAR_H_
