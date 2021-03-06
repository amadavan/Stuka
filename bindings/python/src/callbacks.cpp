//
// Created by Avinash Madavan on 2019-05-09.
//

#include "callbacks.h"

void init_callbacks(py::module &m) {
  py::class_<stuka::util::callback::BaseCallback,
             std::shared_ptr<stuka::util::callback::BaseCallback>>(m, "BaseCallback");

  py::class_<stuka::util::callback::Composite, stuka::util::callback::BaseCallback,
             std::shared_ptr<stuka::util::callback::Composite>>(m, "Composite")
      .def(py::init<const std::vector<std::shared_ptr<stuka::util::callback::BaseCallback>> &>(),
           py::arg("callbacks"));

  py::class_<stuka::util::callback::Function, stuka::util::callback::BaseCallback,
             std::shared_ptr<stuka::util::callback::Function>>(m, "Function")
      .def(py::init<const std::function<void(stuka::OptimizeState)> &>(), py::arg("function"));

  py::class_<stuka::util::callback::Progress, stuka::util::callback::BaseCallback,
             std::shared_ptr<stuka::util::callback::Progress>>(m, "Progress")
      .def(py::init<unsigned int>(), py::arg("n_max_"));

  py::class_<stuka::util::callback::SaveHDF5, stuka::util::callback::BaseCallback,
             std::shared_ptr<stuka::util::callback::SaveHDF5>>(m, "SaveHDF5")
      .def(py::init<const std::string &, const size_t>(), py::arg("filename"), py::arg("n_entries") = 0);
}