//
// Created by Avinash Madavan on 7/16/18.
//

#ifndef STUKA_EXAMPLE_BORRELLI_731_H
#define STUKA_EXAMPLE_BORRELLI_731_H

#include "example.h"

namespace stuka { namespace example {

class Borrelli731 : public DecomposedLinearProgramExample {
 public:
  stuka::dLP::DecomposedLinearProgram gen() override;

  stuka::LP::LinearProgram full() override;

  std::string name() override { return "Borrelli Example 7.31"; }
};
}}

#endif //STUKA_EXAMPLE_BORRELLI_731_H
