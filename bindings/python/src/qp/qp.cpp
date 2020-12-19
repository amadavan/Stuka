//
// Created by Avinash Madavan on 2019-02-05.
//


#include "qp.h"

void init_qp(py::module &m) {

  py::class_<stuka::QP::QuadraticProgram, stuka::QP::PyQuadraticProgram>(m, "QuadraticProgram")
      .def(py::init<Eigen::SparseMatrix<double>, Eigen::VectorXd, Eigen::SparseMatrix<double>, Eigen::VectorXd,
                    Eigen::SparseMatrix<double>, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>(),
           py::arg("Q") = Eigen::SparseMatrix<double>(0, 0),
           py::arg("c") = Eigen::VectorXd(0),
           py::arg("A_ub") = Eigen::SparseMatrix<double>(0, 0),
           py::arg("b_ub") = Eigen::VectorXd(0),
           py::arg("A_eq") = Eigen::SparseMatrix<double>(0, 0),
           py::arg("b_eq") = Eigen::VectorXd(0),
           py::arg("lb") = Eigen::VectorXd(0),
           py::arg("ub") = Eigen::VectorXd(0))
      .def_property("Q",
                    [](stuka::QP::QuadraticProgram &self) {
                      return (self.Q) ? py::cast<Eigen::SparseMatrix<double>>(*self.Q) : py::cast<py::none>(Py_None);
                    }, py::cpp_function([](stuka::QP::QuadraticProgram &self, Eigen::SparseMatrix<double> Q_) {
                                          self.Q = std::make_unique<Eigen::SparseMatrix<double>>(Q_);
                                        },
                                        py::keep_alive<1, 2>()))
      .def_property("c",
                    [](stuka::QP::QuadraticProgram &self) {
                      return (self.c) ? py::cast<Eigen::VectorXd>(*self.c) : py::cast<py::none>(Py_None);
                    },
                    py::cpp_function([](stuka::QP::QuadraticProgram &self,
                                        Eigen::VectorXd c_) { self.c = std::make_unique<Eigen::VectorXd>(c_); },
                                     py::keep_alive<1, 2>()))
      .def_property("A_ub",
                    [](stuka::QP::QuadraticProgram &self) {
                      return (self.A_ub) ? py::cast<Eigen::VectorXd>(*self.A_ub) : py::cast<py::none>(Py_None);
                    },
                    py::cpp_function([](stuka::QP::QuadraticProgram &self,
                                        Eigen::SparseMatrix<double> A_ub_) {
                                       self.A_ub = std::make_unique<Eigen::SparseMatrix<double>>(A_ub_);
                                     },
                                     py::keep_alive<1, 2>()))
      .def_property("b_ub",
                    [](stuka::QP::QuadraticProgram &self) {
                      return (self.b_ub) ? py::cast<Eigen::VectorXd>(*self.b_ub) : py::cast<py::none>(Py_None);
                    },
                    py::cpp_function([](stuka::QP::QuadraticProgram &self,
                                        Eigen::VectorXd b_ub_) {
                                       self.b_ub = std::make_unique<Eigen::VectorXd>(b_ub_);
                                     },
                                     py::keep_alive<1, 2>()))
      .def_property("A_eq",
                    [](stuka::QP::QuadraticProgram &self) {
                      return (self.A_eq) ? py::cast<Eigen::VectorXd>(*self.A_eq) : py::cast<py::none>(Py_None);
                    },
                    py::cpp_function([](stuka::QP::QuadraticProgram &self,
                                        Eigen::SparseMatrix<double> A_eq_) {
                                       self.A_eq = std::make_unique<Eigen::SparseMatrix<double>>(A_eq_);
                                     },
                                     py::keep_alive<1, 2>()))
      .def_property("b_eq",
                    [](stuka::QP::QuadraticProgram &self) {
                      return (self.b_eq) ? py::cast<Eigen::VectorXd>(*self.b_eq) : py::cast<py::none>(Py_None);
                    },
                    py::cpp_function([](stuka::QP::QuadraticProgram &self,
                                        Eigen::VectorXd b_eq_) {
                                       self.b_eq = std::make_unique<Eigen::VectorXd>(b_eq_);
                                     },
                                     py::keep_alive<1, 2>()))
      .def_property("lb",
                    [](stuka::QP::QuadraticProgram &self) {
                      return (self.lb) ? py::cast<Eigen::VectorXd>(*self.lb) : py::cast<py::none>(Py_None);
                    },
                    py::cpp_function([](stuka::QP::QuadraticProgram &self,
                                        Eigen::VectorXd lb_) { self.lb = std::make_unique<Eigen::VectorXd>(lb_); },
                                     py::keep_alive<1, 2>()))
      .def_property("ub",
                    [](stuka::QP::QuadraticProgram &self) {
                      return (self.ub) ? py::cast<Eigen::VectorXd>(*self.ub) : py::cast<py::none>(Py_None);
                    },
                    py::cpp_function([](stuka::QP::QuadraticProgram &self,
                                        Eigen::VectorXd ub_) { self.ub = std::make_unique<Eigen::VectorXd>(ub_); },
                                     py::keep_alive<1, 2>()));
}
