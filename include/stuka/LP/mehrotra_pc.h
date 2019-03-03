//
// Created by Avinash Madavan on 2019-03-02.
//

#ifndef STUKA_LP_MEHROTRA_PC_H
#define STUKA_LP_MEHROTRA_PC_H

#include <Eigen/SparseCholesky>

#include "base_solver.h"
#include "slack_lp.h"

namespace stuka { namespace LP {
  class MehrotraPC : public BaseLPSolver {
  public:
    MehrotraPC(const LinearProgram &lp, const Options &opts);

    void iterate() override;

    bool terminate() override;

    const OptimizeState getState() override;

    BaseLinearProgram &getLP() override;

  private:
    std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solveKKT(Eigen::VectorXd Rc,
                                                                           Eigen::VectorXd Rb,
                                                                           Eigen::VectorXd Rxs,
                                                                           bool resolve = true);

    std::tuple<double, double> computeStepSize(Eigen::VectorXd dx, Eigen::VectorXd ds, double scale = 1.);

    SlackLinearProgram prog_;

    double bc_;
    Eigen::VectorXd e_;

    Eigen::VectorXd x_;
    Eigen::VectorXd y_;
    Eigen::VectorXd s_;

    double mu_;
    Eigen::VectorXd Rc_;
    Eigen::VectorXd Rb_;
    Eigen::VectorXd Rxs_;
    double error_;
    const double eps_;

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver_;
  };
}}

#endif //STUKA_LP_MEHROTRA_PC_H
