//
// Created by Avinash Madavan on 1/15/19.
//

#ifndef STUKA_dLP_CRITICAL_REGION_H
#define STUKA_dLP_CRITICAL_REGION_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace stuka { namespace dLP {
struct CriticalRegion {
  // Cost structure (alpha * x + beta)
  Eigen::VectorXd alpha;
  double beta = 0;

  // Boundaries (A * x + A_add * x_add <= b)
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;

  size_t n_add = 0;
//  Eigen::VectorXd c_add;
  Eigen::SparseMatrix<double> A_add;

  /* Default constructor
   *
   */
  CriticalRegion() = default;

  /* Copy constructor
   *
   */
  CriticalRegion(const CriticalRegion &cr) = default;

  /* Determine whether point lies within critical region
   *
   * Returns whether the provided point lies within the critical region by evaluating whether the
   * inequality is satisfied. Assumes that n_add = 0.
   *
   * TODO: incorporate potential dual degeneracy (n_add > 0)
   */
  bool in(const Eigen::VectorXd &x) const;
};
}}

#endif //STUKA_dLP_CRITICAL_REGION_H
