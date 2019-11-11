//
// Created by Avinash Madavan on 2019-09-16.
//

#ifndef STUKA_STOCHASTIC_MEASURE_EXPECTED_VALUE_H_
#define STUKA_STOCHASTIC_MEASURE_EXPECTED_VALUE_H_

#include "../program.h"
#include "measure.h"

namespace stuka { namespace stochastic {
class ExpectedValue : public ProgramMeasure {
 public:
  ExpectedValue() {}
  ~ExpectedValue() {}

  double f(const Program &p, const Eigen::VectorXd &x, const Eigen::VectorXd &xi) override {
    return p.f(x, xi);
  }
  Eigen::VectorXd g(const Program &p, const Eigen::VectorXd &x, const Eigen::VectorXd &xi) override {
    return p.g(x, xi);
  }
  Eigen::VectorXd df(const Program &p, const Eigen::VectorXd &x, const Eigen::VectorXd &xi) override {
    return p.df(x, xi);
  }
  Eigen::MatrixXd dg(const Program &p, const Eigen::VectorXd &x, const Eigen::VectorXd &xi) override {
    return p.dg(x, xi);
  }
  Eigen::VectorXd projX(const Program &p, const Eigen::VectorXd &x) override {
    return p.projX(x);
  }
};
}}

#endif //STUKA_STOCHASTIC_MEASURE_EXPECTED_VALUE_H_
