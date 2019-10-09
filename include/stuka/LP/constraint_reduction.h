//
// Created by Avinash Madavan on 10/8/19.
//

#ifndef STUKA_UTIL_CONSTRAINT_REDUCTION_H_
#define STUKA_UTIL_CONSTRAINT_REDUCTION_H_

#include <stddef.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "lp.h"

namespace stuka { namespace LP {

class ConstraintReduction {
 public:
  static void bounds(const LinearProgram &prog);
};

}}

#endif //STUKA_UTIL_CONSTRAINT_REDUCTION_H_
