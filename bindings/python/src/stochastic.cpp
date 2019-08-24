//
// Created by Avinash Madavan on 2019-05-07.
//

#include "stochastic.h"

void init_stochastic(py::module &m) {
  py::class_<stuka::stochastic::Program>(m, "StochasticProgram")
      .def(py::init<>())
      .def_readwrite("f", &stuka::stochastic::Program::f)
      .def_readwrite("g", &stuka::stochastic::Program::g)
      .def_readwrite("df", &stuka::stochastic::Program::df)
      .def_readwrite("dg", &stuka::stochastic::Program::dg)
      .def_readwrite("projX", &stuka::stochastic::Program::projX)
      .def_readwrite("sample", &stuka::stochastic::Program::sample)
      .def_readwrite("alpha", &stuka::stochastic::Program::alpha)
      .def_readwrite("beta", &stuka::stochastic::Program::beta)
      .def_readwrite("gmax", &stuka::stochastic::Program::gmax);
}