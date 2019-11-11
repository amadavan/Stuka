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
      .def_readwrite("sample", &stuka::stochastic::Program::sample);

  py::module measure = m.def_submodule("measure");
  py::class_ < stuka::stochastic::ExpectedValue,
      std::shared_ptr < stuka::stochastic::ExpectedValue >> (measure, "ExpectedValue")
          .def(py::init<>());
  py::class_ < stuka::stochastic::ConditionalValueAtRisk,
      std::shared_ptr < stuka::stochastic::ConditionalValueAtRisk >> (measure, "ConditionalValueAtRisk")
          .def(py::init<size_t, double, Eigen::VectorXd, Eigen::VectorXd>(),
               py::arg("n_var"),
               py::arg("alpha") = 0.,
               py::arg("beta") = Eigen::VectorXd(),
               py::arg("g_max"));
}