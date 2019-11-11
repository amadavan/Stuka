//
// Created by Avinash Madavan on 2019-08-26.
//

#ifndef STUKA_STOCHASTIC_MEASURE_MEASURE_H_
#define STUKA_STOCHASTIC_MEASURE_MEASURE_H_

#include "../program.h"

namespace stuka { namespace stochastic {

/* Apply a measure to a stochastic program.
 *
 * Various measures can be applied to stochastic programs, including expectation, CVaR, EVaR, or other risk measures.
 * Each of these measures augments the original stochastic program, which we handle by introducing the notation of
 * measure and any associated additional variables that are required.
 */
class ProgramMeasure {
 public:
  ProgramMeasure() {}
  virtual ~ProgramMeasure() = default;

  virtual double f(const Program &, const Eigen::VectorXd &, const Eigen::VectorXd &) = 0;
  virtual Eigen::VectorXd g(const Program &, const Eigen::VectorXd &, const Eigen::VectorXd &) = 0;
  virtual Eigen::VectorXd df(const Program &, const Eigen::VectorXd &, const Eigen::VectorXd &) = 0;
  virtual Eigen::MatrixXd dg(const Program &, const Eigen::VectorXd &, const Eigen::VectorXd &) = 0;
  virtual Eigen::VectorXd projX(const Program &, const Eigen::VectorXd &) = 0;
  virtual Eigen::VectorXd augmentState(Eigen::VectorXd &x) { return x; }
};

}}

#endif //STUKA_STOCHASTIC_MEASURE_MEASURE_H_
