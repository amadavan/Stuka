//
// Created by Avinash Madavan on 1/21/19.
//

#ifndef STUKA_PYTHON_LP_H
#define STUKA_PYTHON_LP_H

#include <stuka/LP/lp.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace stuka { namespace LP {
  class PyLinearProgram : public LinearProgram {
  public:
    using LinearProgram::LinearProgram;

    PyLinearProgram(Eigen::VectorXd _c,
                    Eigen::SparseMatrix<double> _A_ub,
                    Eigen::VectorXd _b_ub,
                    Eigen::SparseMatrix<double> _A_eq,
                    Eigen::VectorXd _b_eq,
                    Eigen::VectorXd _lb,
                    Eigen::VectorXd _ub) {
      if (_c.size() > 0) c = std::make_unique<Eigen::VectorXd>(_c);
      if(_b_ub.size() > 0) {
        A_ub = std::make_unique<Eigen::SparseMatrix<double>>(_A_ub);
        b_ub = std::make_unique<Eigen::VectorXd>(_b_ub);
      }
      if (_b_eq.size() > 0) {
        A_eq = std::make_unique<Eigen::SparseMatrix<double>>(_A_eq);
        b_eq = std::make_unique<Eigen::VectorXd>(_b_eq);
      }
      if (_lb.size() > 0) lb = std::make_unique<Eigen::VectorXd>(_lb);
      if (_ub.size() > 0) ub = std::make_unique<Eigen::VectorXd>(_ub);
    }
  };
}}

namespace py = pybind11;

void init_lp(py::module &);

#endif //STUKA_PYTHON_LP_H
