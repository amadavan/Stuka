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
    const Program &prog_;   // Main program data

    Eigen::VectorXd x_;     // Primal variable
    Eigen::VectorXd z_;     // Dual variable with inequality constraint
    Eigen::VectorXd xbar_;  // Average primal variable
    Eigen::VectorXd zbar_;  // Average dual variable with inequality constraint

    size_t n_;              // Number of variables
    size_t m_;              // Number of constraints
    size_t m_g_;            // Number of stochastic constraints
    size_t m_h_;            // Number of deterministic constraints

    double step_;           // Step size

    double alphap_;         // 1/(1 - alpha)
    Eigen::VectorXd betap_; // 1/(1 - beta)

    size_t nit_;            // iteration

    std::function<double(Eigen::VectorXd &, Eigen::VectorXd &)> f_;
    std::function<Eigen::VectorXd(Eigen::VectorXd &, Eigen::VectorXd &)> g_;
    std::function<Eigen::VectorXd(Eigen::VectorXd &, Eigen::VectorXd &)> df_;
    std::function<Eigen::MatrixXd(Eigen::VectorXd &, Eigen::VectorXd &)> dg_;
    std::function<Eigen::VectorXd(Eigen::VectorXd &)> projX_;

    double f_cvar(Eigen::VectorXd &xu, Eigen::VectorXd &xi);

    Eigen::VectorXd g_cvar(Eigen::VectorXd &xu, Eigen::VectorXd &xi);

    Eigen::VectorXd df_cvar(Eigen::VectorXd &xu, Eigen::VectorXd &xi);

    Eigen::MatrixXd dg_cvar(Eigen::VectorXd &xu, Eigen::VectorXd &xi);

    Eigen::VectorXd projX_cvar(Eigen::VectorXd &xu);
  };
}}

#endif //STUKA_STOCHASTIC_EXPECTED_PD_H
