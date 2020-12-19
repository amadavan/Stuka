//
// Created by Avinash Madavan on 1/10/19.
//

#ifndef STUKA_dLP_DECOMPOSED_LP_H
#define STUKA_dLP_DECOMPOSED_LP_H

#include <vector>
#include <memory>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../constants.h"
#include "../util/dense_ops.h"
#include "../util/sparse_ops.h"

namespace stuka { namespace dLP {
struct DecomposedLinearProgram {
  std::vector<std::unique_ptr<Eigen::VectorXd>> c;
  std::vector<std::unique_ptr<Eigen::SparseMatrix<double>>> A_ub;
  std::vector<std::unique_ptr<Eigen::VectorXd>> b_ub;
  std::vector<std::unique_ptr<Eigen::SparseMatrix<double>>> C_ub;
  std::vector<std::unique_ptr<Eigen::SparseMatrix<double>>> A_eq;
  std::vector<std::unique_ptr<Eigen::VectorXd>> b_eq;
  std::vector<std::unique_ptr<Eigen::SparseMatrix<double>>> C_eq;
  std::vector<std::unique_ptr<Eigen::VectorXd>> lb;
  std::vector<std::unique_ptr<Eigen::VectorXd>> ub;

  DecomposedLinearProgram() = default;

  explicit DecomposedLinearProgram(size_t n_sub) {
    c = std::vector<std::unique_ptr<Eigen::VectorXd>>(n_sub + 1);
    A_ub = std::vector<std::unique_ptr<Eigen::SparseMatrix<double>>>(n_sub + 1);
    b_ub = std::vector<std::unique_ptr<Eigen::VectorXd>>(n_sub + 1);
    C_ub = std::vector<std::unique_ptr<Eigen::SparseMatrix<double>>>(n_sub);
    A_eq = std::vector<std::unique_ptr<Eigen::SparseMatrix<double>>>(n_sub + 1);
    b_eq = std::vector<std::unique_ptr<Eigen::VectorXd>>(n_sub + 1);
    C_eq = std::vector<std::unique_ptr<Eigen::SparseMatrix<double>>>(n_sub);
    lb = std::vector<std::unique_ptr<Eigen::VectorXd>>(n_sub + 1);
    ub = std::vector<std::unique_ptr<Eigen::VectorXd>>(n_sub + 1);
  }

  DecomposedLinearProgram(
      std::vector<std::unique_ptr<Eigen::VectorXd>> _c,
      std::vector<std::unique_ptr<Eigen::SparseMatrix<double>>> _A_ub,
      std::vector<std::unique_ptr<Eigen::VectorXd>> _b_ub,
      std::vector<std::unique_ptr<Eigen::SparseMatrix<double>>> _C_ub,
      std::vector<std::unique_ptr<Eigen::SparseMatrix<double>>> _A_eq,
      std::vector<std::unique_ptr<Eigen::VectorXd>> _b_eq,
      std::vector<std::unique_ptr<Eigen::SparseMatrix<double>>> _C_eq,
      std::vector<std::unique_ptr<Eigen::VectorXd>> _lb,
      std::vector<std::unique_ptr<Eigen::VectorXd>> _ub) : DecomposedLinearProgram(_c.size() - 1) {
    for (size_t i = 0; i < _c.size(); ++i) {
      c[i] = util::DenseOps::unique_copy(_c[i]);
      A_ub[i] = util::SparseOps::unique_copy(_A_ub[i]);
      b_ub[i] = util::DenseOps::unique_copy(_b_ub[i]);
      A_eq[i] = util::SparseOps::unique_copy(_A_eq[i]);
      b_eq[i] = util::DenseOps::unique_copy(_b_eq[i]);
      lb[i] = util::DenseOps::unique_copy(_lb[i]);
      ub[i] = util::DenseOps::unique_copy(_ub[i]);
      if (i < _c.size() - 1) {
        C_ub[i] = util::SparseOps::unique_copy(_C_ub[i]);
        C_eq[i] = util::SparseOps::unique_copy(_C_eq[i]);
      }
    }
  }

  void reduceConstraints(ConstraintReductionMethods method = BOUNDS) {
    switch (method) {
      case BOUNDS: ConstraintReductionBounds();
        break;
    }
  };

 private:
  // Hack for cleaner python implementation
  virtual std::string name() { return "DecomposedLinearProgram"; }

  void ConstraintReductionBounds();
};
}}

#endif //STUKA_dLP_DECOMPOSED_LP_H
