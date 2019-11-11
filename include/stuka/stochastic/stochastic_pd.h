//
// Created by Avinash Madavan on 2019-05-07.
//

#ifndef STUKA_STOCHASTIC_EXPECTED_PD_H
#define STUKA_STOCHASTIC_EXPECTED_PD_H

#include <algorithm>
#include <functional>

#include "base_solver.h"

namespace stuka { namespace stochastic {
class StochasticPrimalDual : public BaseStochasticSolver {
 public:
  StochasticPrimalDual(const Program &prog, const Options &opts);

  void iterate() override;

  bool terminate() override;

  const OptimizeState getState() override;

 private:
  const Program &prog_;                             // Main program data
  const std::shared_ptr<ProgramMeasure> measure_;   // Measure on which to optimize

  Eigen::VectorXd x_;     // Primal variable
  Eigen::VectorXd z_;     // Dual variable with inequality constraint
  Eigen::VectorXd x_bar_;  // Average primal variable
  Eigen::VectorXd z_bar_;  // Average dual variable with inequality constraint

  size_t n_;      // Number of variables
  size_t m_;      // Number of constraints
  size_t m_g_;    // Number of stochastic constraints
  size_t m_h_;    // Number of deterministic constraints

  double step_;     // Step size
  size_t nit_;      // iteration
};
}}

#endif //STUKA_STOCHASTIC_EXPECTED_PD_H
