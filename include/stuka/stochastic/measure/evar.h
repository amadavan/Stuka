//
// Created by Avinash Madavan on 10/8/19.
//

#ifndef STUKA_STOCHASTIC_MEASURE_EVAR_H_
#define STUKA_STOCHASTIC_MEASURE_EVAR_H_

#include <math.h>

#include "g_entropic.h"

namespace stuka { namespace stochastic {

class EntropicConjugateFunction : public gConjugateFunction {
 public:
  double g(const double &x) override {
    return exp(static_cast<long double>(x - 1));
  }
  double dg(const double &x) override {
    return exp(static_cast<long double>(x - 1));
  }
};

using EntropicValueAtRisk = gEntropicRisk<EntropicConjugateFunction>;

}}

#endif //STUKA_STOCHASTIC_MEASURE_EVAR_H_
