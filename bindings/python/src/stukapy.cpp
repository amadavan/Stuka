//
// Created by Avinash Madavan on 1/21/19.
//

#include <string>
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <stuka/stuka.h>

#include "dlp.h"
#include "lp.h"
#include "qp.h"

namespace py = pybind11;

namespace pybind11 { namespace detail {
  template<>
  struct type_caster<stuka::OptimizeState> {
  public:
  PYBIND11_TYPE_CASTER(stuka::OptimizeState, _("OptimizeState"));

    bool load(handle src, bool) {
      return false;
    }

    static handle cast(stuka::OptimizeState res, return_value_policy /* policy */, handle /* parent */) {
      py::dict ret = py::dict(
          "x"_a = res.x,
          "dual_ub"_a = res.dual_ub,
          "dual_eq"_a = res.dual_eq,
          "fun"_a = res.fun,
          "error"_a = res.error,
          "status"_a = res.status,
          "nit"_a = res.nit,
          "nit_sub"_a = res.nit_sub,
          "runtime"_a = res.runtime
      );

      return ret.release();
    }
  };

}}

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

  m.attr("inf") = stuka::INF;

  py::enum_<stuka::Solver>(m, "solver")
#ifdef BUILD_GUROBI
      .value("GUROBI", stuka::Solver::GUROBI)
#endif
      .value("BENDER", stuka::Solver::BENDER)
      .value("CRE", stuka::Solver::CRE)
      .value("MPC", stuka::Solver::MPC);

  py::class_<stuka::Options>(m, "Options")
      .def(py::init<>())
      .def_readwrite("x0", &stuka::Options::x0)
      .def_readwrite("tol", &stuka::Options::tol)
      .def_readwrite("max_iter", &stuka::Options::max_iter)
      .def_readwrite("dual_ub0", &stuka::Options::dual_ub0)
      .def_readwrite("dual_eq0", &stuka::Options::dual_eq0)
      .def_readwrite("lp_solver", &stuka::Options::lp_solver)
      .def_readwrite("qp_solver", &stuka::Options::qp_solver)
      .def_readwrite("dlp_solver", &stuka::Options::dlp_solver);

  init_lp(m);
  init_qp(m);

  py::class_<stuka::dLP::DecomposedLinearProgram, stuka::dLP::PyDecomposedLinearProgram>(m, "DecomposedLinearProgram")
      .def(
          py::init<std::vector<Eigen::VectorXd>, std::vector<Eigen::SparseMatrix<double>>, std::vector<Eigen::VectorXd>,
              std::vector<Eigen::SparseMatrix<double>>, std::vector<Eigen::SparseMatrix<double>>,
              std::vector<Eigen::VectorXd>, std::vector<Eigen::SparseMatrix<double>>, std::vector<Eigen::VectorXd>,
              std::vector<Eigen::VectorXd>>(),
          py::arg("c"), py::arg("A_ub"), py::arg("b_ub"), py::arg("C_ub"), py::arg("A_eq"), py::arg("b_eq"),
          py::arg("C_eq"), py::arg("lb"), py::arg("ub"));

  m.def("linprog", py::overload_cast<const stuka::LP::LinearProgram &, const stuka::Options &>(&stuka::util::linprog),
        R"pbdoc( solve a linear program )pbdoc", py::arg("LinearProgram"), py::arg("Options") = stuka::Options());

  m.def("quadprog", &stuka::util::quadprog, R"pbdoc( solve a quadratic program )pbdoc", py::arg("QuadraticProgram"),
        py::arg("Options") = stuka::Options());

  m.def("linprog",
        py::overload_cast<const stuka::dLP::DecomposedLinearProgram &, const stuka::Options &>(&stuka::util::linprog),
        R"pbdoc( solve a decomposed linear program )pbdoc",
        py::arg("DecomposedLinearProgram"), py::arg("Options") = stuka::Options());

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}