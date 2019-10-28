//
// Created by Avinash Madavan on 2019-05-07.
//

#ifndef STUKA_STOCHASTIC_PROGRAM_H
#define STUKA_STOCHASTIC_PROGRAM_H

#include <Eigen/Core>
#include <functional>

namespace stuka { namespace stochastic {
/* A General Class of problems
 *
 * Defined the general optimization problem:
 *    min  f(x),
 *    s.t. g(x) <= 0
 */
struct Program {
  std::function<double(const Eigen::VectorXd &,
                       const Eigen::VectorXd &)> f;                     // Function value given iterate/sample
  std::function<Eigen::VectorXd(const Eigen::VectorXd &,
                                const Eigen::VectorXd &)> g;            // Constraint values given iterate/sample
  std::function<Eigen::VectorXd(const Eigen::VectorXd &)> h;            // Deterministic constraint values
  std::function<Eigen::VectorXd(const Eigen::VectorXd &,
                                const Eigen::VectorXd &)> df;           // Function gradient given iterate/sample
  std::function<Eigen::MatrixXd(const Eigen::VectorXd &,
                                const Eigen::VectorXd &)> dg;           // Stochastic Constraint gradients
  std::function<Eigen::MatrixXd(const Eigen::VectorXd &)> dh;           // Deterministic constraint gradients
  std::function<Eigen::VectorXd(const Eigen::VectorXd &)> projX;        // Projection of iterate onto feasible set

  std::function<Eigen::VectorXd()> sample;                              // Acquire a random sample

  double alpha = 0.;                                                    // Function CVaR parameter
  Eigen::VectorXd beta = Eigen::VectorXd(0);                            // Constraint CVaR parameter
  Eigen::VectorXd gmax = Eigen::VectorXd(0);                            // Max constraint value for PDSS

  // Hack for cleaner python implementation
  virtual std::string name() { return "StochasticProgram"; }
};
}}

#endif //STUKA_PROGRAM_H
