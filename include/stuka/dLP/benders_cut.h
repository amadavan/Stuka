//
// Created by Avinash Madavan on 1/15/19.
//

#ifndef STUKA_dLP_BENDERS_CUT_H
#define STUKA_dLP_BENDERS_CUT_H

#include <Eigen/Core>

namespace stuka { namespace dLP {
  struct BendersCut {
    Eigen::VectorXd a;
    double b;
  };
}}

#endif //STUKA_dLP_BENDERS_CUT_H
