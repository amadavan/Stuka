//
// Created by Avinash Madavan on 2019-01-29.
//

#ifndef STUKA_PYTHON_QP_H
#define STUKA_PYTHON_QP_H

#include <stuka/QP/qp.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace stuka { namespace QP {
  class PyQuadraticProgram : public QuadraticProgram {
  public:
    using QuadraticProgram::QuadraticProgram;

    PyQuadraticProgram(Eigen::SparseMatrix<double> _Q,
                       Eigen::VectorXd _c,
                       Eigen::SparseMatrix<double> _A_ub,
                       Eigen::VectorXd _b_ub,
                       Eigen::SparseMatrix<double> _A_eq,
                       Eigen::VectorXd _b_eq,
                       Eigen::VectorXd _lb,
                       Eigen::VectorXd _ub) {
      if (_Q.rows() > 0 && _Q.cols() > 0) Q = std::make_shared<Eigen::SparseMatrix<double>>(_Q);
      if (_c.size() > 0) c = std::make_shared<Eigen::VectorXd>(_c);
      if (_b_ub.size() > 0) {
        A_ub = std::make_shared<Eigen::SparseMatrix<double>>(_A_ub);
        b_ub = std::make_shared<Eigen::VectorXd>(_b_ub);
      }
      if (_b_eq.size() > 0) {
        A_eq = std::make_shared<Eigen::SparseMatrix<double>>(_A_eq);
        b_eq = std::make_shared<Eigen::VectorXd>(_b_eq);
      }
      if (_lb.size() > 0) lb = std::make_shared<Eigen::VectorXd>(_lb);
      if (_ub.size() > 0) ub = std::make_shared<Eigen::VectorXd>(_ub);
    }
  };
}}

namespace py = pybind11;

void init_qp(py::module &);

#endif //STUKA_PYTHON_QP_H
