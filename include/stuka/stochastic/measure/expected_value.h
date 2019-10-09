//
// Created by Avinash Madavan on 2019-09-16.
//

#ifndef STUKA_STOCHASTIC_MEASURE_EXPECTED_VALUE_H_
#define STUKA_STOCHASTIC_MEASURE_EXPECTED_VALUE_H_

#include "measure.h"

namespace stuka { namespace stochastic {
  class ExpectedValue : public MeasurableProgram {
   public:
    ExpectedValue(const Program &prog) : MeasurableProgram(prog) {}
  };
}}

#endif //STUKA_STOCHASTIC_MEASURE_EXPECTED_VALUE_H_
