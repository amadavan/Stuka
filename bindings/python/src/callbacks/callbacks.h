//
// Created by Avinash Madavan on 2019-05-09.
//

#ifndef STUKA_PYTHON_CALLBACKS_H
#define STUKA_PYTHON_CALLBACKS_H

#include <stuka/optimize_state.h>
#include <stuka/util/callback/base_callback.h>
#include <stuka/util/callback/function.h>
#include <stuka/util/callback/save_hdf5.h>
#include <stuka/util/callback/composite.h>
#include <stuka/util/callback/progress.h>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

namespace py = pybind11;

void init_callbacks(py::module &);

#endif //STUKA_PYTHON_CALLBACKS_H
