//
// Created by Avinash Madavan on 2019-05-08.
//

#include <stuka/util/callback/function.h>

stuka::util::callback::Function::Function(const std::function<void(OptimizeState)> func) : func_(func) {};

void stuka::util::callback::Function::callback(stuka::OptimizeState state) {
  return func_(state);
}
