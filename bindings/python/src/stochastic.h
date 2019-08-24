//
// Created by Avinash Madavan on 2019-05-07.
//

#ifndef STUKA_PYTHON_STOCHASTIC_H
#define STUKA_PYTHON_STOCHASTIC_H

#include <stuka/stochastic/program.h>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_stochastic(py::module &);

#endif //STUKA_PYTHON_STOCHASTIC_H
