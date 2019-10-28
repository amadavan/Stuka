//
// Created by Avinash Madavan on 7/16/18.
//

#ifndef STUKA_EXAMPLE_EXAMPLE_H
#define STUKA_EXAMPLE_EXAMPLE_H

#include <string>

#include <stuka/constants.h>
#include <stuka/LP/lp.h>
#include <stuka/QP/qp.h>
#include <stuka/dLP/decomposed_lp.h>

namespace stuka { namespace example {

class LinearProgramExample {
 public:
  virtual ~LinearProgramExample() = default;

  virtual stuka::LP::LinearProgram gen() = 0;

  virtual std::string name() = 0;
};

class QuadraticProgramExample {
 public:
  virtual ~QuadraticProgramExample() = default;

  virtual stuka::QP::QuadraticProgram gen() = 0;

  virtual std::string name() = 0;
};

class DecomposedLinearProgramExample {
 public:
  virtual ~DecomposedLinearProgramExample() = default;

  virtual stuka::dLP::DecomposedLinearProgram gen() = 0;

  virtual stuka::LP::LinearProgram full() = 0;

  virtual std::string name() = 0;
};

}}

#endif //STUKA_EXAMPLE_EXAMPLE_H
