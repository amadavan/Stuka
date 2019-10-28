//
// Created by Avinash Madavan on 12/23/18.
//

#ifndef STUKA_LP_LP_H
#define STUKA_LP_LP_H

#include <memory>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../constants.h"

namespace stuka { namespace LP {
struct LinearProgram {
  std::shared_ptr<Eigen::VectorXd> c;
  std::shared_ptr<Eigen::SparseMatrix<double>> A_ub;
  std::shared_ptr<Eigen::VectorXd> b_ub;
  std::shared_ptr<Eigen::SparseMatrix<double>> A_eq;
  std::shared_ptr<Eigen::VectorXd> b_eq;
  std::shared_ptr<Eigen::VectorXd> lb;
  std::shared_ptr<Eigen::VectorXd> ub;

  void reduceConstraints(ConstraintReductionMethods method = BOUNDS) {
    switch (method) {
      case BOUNDS: ConstraintReductionBounds();
        break;
    }
  };

  // Hack for cleaner python implementation
  virtual std::string name() { return "LinearProgram"; }

 private:
  void ConstraintReductionBounds();
};
}}

#endif //STUKA_LP_LP_H
