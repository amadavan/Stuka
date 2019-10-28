//
// Created by Avinash Madavan on 7/18/18.
//

#ifndef SWITCHBACK_EXAMPLE_SIMPLE_QP_H
#define SWITCHBACK_EXAMPLE_SIMPLE_QP_H

#include "example.h"

namespace stuka { namespace example {

// Example from https://www.mathworks.com/help/optim/ug/quadprog.html
class SimpleQP : public QuadraticProgramExample {
 public:
  QP::QuadraticProgram gen() override;

  std::string name() override;
};

}}

#endif //SWITCHBACK_EXAMPLE_SIMPLE_QP_H
