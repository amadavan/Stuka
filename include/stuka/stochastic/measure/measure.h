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
class MeasurableProgram : public Program {
 public:
  MeasurableProgram(const Program &prog) {
    f = prog.f;
    g = prog.g;
    h = prog.h;
    df = prog.df;
    dg = prog.dg;
    dh = prog.dh;
    projX = prog.projX;
    sample = prog.sample;
  }
  virtual ~MeasurableProgram() {}
};

}}

#endif //STUKA_STOCHASTIC_MEASURE_MEASURE_H_
