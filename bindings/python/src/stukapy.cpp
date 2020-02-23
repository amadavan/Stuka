//
// Created by Avinash Madavan on 1/21/19.
//

#include <string>
#include <memory>
#include <variant>
#include <iostream>
#include <iomanip>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <stuka/stuka.h>

#include "lp/lp.h"
#include "qp/qp.h"
#include "dlp/dlp.h"
#include "stochastic/stochastic.h"
#include "callbacks/callbacks.h"

namespace py = pybind11;

std::ostream &operator<<(std::ostream &os, const stuka::OptimizeState &state) {
  os << std::right << std::setw(8) << std::setfill(' ') << "x" << ": " << std::left << state.x.transpose() << std::endl;
  os << std::right << std::setw(8) << std::setfill(' ') << "dual_ub" << ": " << std::left << state.dual_ub.transpose()
     << std::endl;
  os << std::right << std::setw(8) << std::setfill(' ') << "dual_eq" << ": " << std::left << state.dual_eq.transpose()
     << std::endl;
  os << std::right << std::setw(8) << std::setfill(' ') << "fun" << ": " << std::left << state.fun << std::endl;
  os << std::right << std::setw(8) << std::setfill(' ') << "error" << ": " << std::left << state.error << std::endl;
  os << std::right << std::setw(8) << std::setfill(' ') << "status" << ": " << std::left << state.status << std::endl;
  os << std::right << std::setw(8) << std::setfill(' ') << "nit" << ": " << std::left << state.nit << std::endl;
  os << std::right << std::setw(8) << std::setfill(' ') << "nit_sub" << ": " << std::left << state.nit_sub << std::endl;
  os << std::right << std::setw(8) << std::setfill(' ') << "runtime" << ": " << std::left << state.runtime << std::endl;
  return os;
}

PYBIND11_MODULE(stukapy, m) {
  m.doc() = R"pbdoc(
            Example Module
            --------------
            .. current module:: pymnopt
            .. autosummary::
                :toctree: _generate
        )pbdoc";

//  py::module mlp = m.def_submodule("lp", "linear programming submodule");
//  py::module mqp = m.def_submodule("qp", "quadratic programming submodule");
  py::module mcallback = m.def_submodule("callback", "callbacks submodule");

  m.attr("inf") = stuka::INF;

  py::enum_<stuka::Solver>(m, "solver")
#ifdef ENABLE_GUROBI
      .value("GUROBI", stuka::Solver::GUROBI)
#endif
      .value("BENDER", stuka::Solver::BENDER)
      .value("CRE", stuka::Solver::CRE)
      .value("MPC", stuka::Solver::MPC)
      .value("PDSS", stuka::Solver::PDSS)
      .value("PDSS2", stuka::Solver::PDSS2);

  py::enum_<stuka::ConstraintReductionMethods>(m, "constraintReductionMethods")
      .value("BOUNDS", stuka::ConstraintReductionMethods::BOUNDS);

  py::class_<stuka::OptimizeState>(m, "OptimizeState")
      .def(py::init<>())
      .def_readonly("x", &stuka::OptimizeState::x)
      .def_readonly("dual_ub", &stuka::OptimizeState::dual_ub)
      .def_readonly("dual_eq", &stuka::OptimizeState::dual_eq)
      .def_readonly("fun", &stuka::OptimizeState::fun)
      .def_readonly("error", &stuka::OptimizeState::error)
      .def_readonly("status", &stuka::OptimizeState::status)
      .def_readonly("nit", &stuka::OptimizeState::nit)
      .def_readonly("nit_sub", &stuka::OptimizeState::nit_sub)
      .def_readonly("runtime", &stuka::OptimizeState::runtime)
      .def("__repr__", [](const stuka::OptimizeState &state) {
        std::stringstream ss;
        ss << state << std::endl;
        return ss.str();
      })
      .def("__getitem__", [](const stuka::OptimizeState &state, const std::string &key) {
             std::variant<Eigen::VectorXd, double, size_t, int, std::string> ret = "Key not found.";;
             if (key == "x") ret = state.x;
             if (key == "dual_ub") ret = state.dual_ub;
             if (key == "dual_eq") ret = state.dual_eq;
             if (key == "fun") ret = state.fun;
             if (key == "error") ret = state.error;
             if (key == "status") ret = state.status;
             if (key == "nit") ret = state.nit;
             if (key == "nit_sub") ret = state.nit_sub;
             if (key == "runtime") ret = state.runtime;
             return ret;
           },
           py::is_operator()
      );

  py::class_<stuka::Options>(m, "Options")
      .def(py::init<>())
      .def_readwrite("x0", &stuka::Options::x0)
      .def_readwrite("tol", &stuka::Options::tol)
      .def_readwrite("max_iter", &stuka::Options::max_iter)
      .def_readwrite("step", &stuka::Options::step)
      .def_readwrite("dual_ub0", &stuka::Options::dual_ub0)
      .def_readwrite("dual_eq0", &stuka::Options::dual_eq0)
      .def_readwrite("measure", &stuka::Options::measure)
      .def_readwrite("lp_solver", &stuka::Options::lp_solver)
      .def_readwrite("qp_solver", &stuka::Options::qp_solver)
      .def_readwrite("dlp_solver", &stuka::Options::dlp_solver)
.def_readwrite("cre_step", &stuka::Options::cre_step)
      .def_readwrite("callback", &stuka::Options::callback);

  init_lp(m);
  init_qp(m);
  init_dlp(m);
  init_stochastic(m);
  init_callbacks(mcallback);

  m.def("linprog", py::overload_cast<const stuka::LP::LinearProgram &, const stuka::Options &>(&stuka::util::linprog),
        R"pbdoc( solve a linear program )pbdoc",
        py::arg("LinearProgram"),
        py::arg("Options") = stuka::Options()
  );

  m.def("quadprog", &stuka::util::quadprog, R"pbdoc( solve a quadratic program )pbdoc",
        py::arg("QuadraticProgram"),
        py::arg("Options") = stuka::Options());

  m.def("linprog",
        py::overload_cast<const stuka::dLP::DecomposedLinearProgram &, const stuka::Options &>(&stuka::util::linprog),
        R"pbdoc( solve a decomposed linear program )pbdoc",
        py::arg("DecomposedLinearProgram"),
        py::arg("Options") = stuka::Options()
  );

  m.def("stochastic", &stuka::util::stochastic, R"pbdoc( solve a CVaR-stochastic convex problem )pbdoc",
        py::arg("Program"),
        py::arg("Options") = stuka::Options()
  );

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}