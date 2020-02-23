//
// Created by Avinash Madavan on 1/15/19.
//

#include <stuka/dLP/critical_region.h>

bool stuka::dLP::CriticalRegion::in(const Eigen::VectorXd &x) const {
  return (A.cols() == x.size()) && ((b - A * x).array() >= -1e-8).all();
}
