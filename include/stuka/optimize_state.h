//
// Created by Avinash Madavan on 1/10/19.
//

#ifndef STUKA_OPTIMIZE_STATE_H
#define STUKA_OPTIMIZE_STATE_H

#include <Eigen/Core>

namespace stuka {

struct OptimizeState {
  Eigen::VectorXd x;
  Eigen::VectorXd dual_ub;
  Eigen::VectorXd dual_eq;
  Eigen::VectorXd dual_x_lb;
  Eigen::VectorXd dual_x_ub;

  double fun = 0;
  double error = 0;
  int status = 0;
  size_t nit = 0;
  size_t nit_sub = 0;

  double runtime = 0;
};

}

#endif //STUKA_OPTIMIZE_STATE_H
