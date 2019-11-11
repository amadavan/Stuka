//
// Created by Avinash Madavan on 2019-05-07.
//

#include <stuka/constants.h>
#include <stuka/stochastic/measure/expected_value.h>
#include "dlp.h"

void init_dlp(py::module &m) {

  py::class_<stuka::dLP::DecomposedLinearProgram, stuka::dLP::PyDecomposedLinearProgram>(m, "DecomposedLinearProgram")
      .def(py::init<std::vector<Eigen::VectorXd>,
                    std::vector<Eigen::SparseMatrix<double >>,
                    std::vector<Eigen::VectorXd>,
                    std::vector<Eigen::SparseMatrix<double >>,
                    std::vector<Eigen::SparseMatrix<double >>,
                    std::vector<Eigen::VectorXd>,
                    std::vector<Eigen::SparseMatrix<double >>,
                    std::vector<Eigen::VectorXd>,
                    std::vector<Eigen::VectorXd >>(),
           py::arg("c"), py::arg("A_ub"), py::arg("b_ub"),
           py::arg("C_ub"), py::arg("A_eq"), py::arg("b_eq"),
           py::arg("C_eq"), py::arg("lb"), py::arg("ub"))
      .def("reduceConstraints",
           &stuka::dLP::DecomposedLinearProgram::reduceConstraints,
           py::arg("method") = stuka::ConstraintReductionMethods::BOUNDS);

}