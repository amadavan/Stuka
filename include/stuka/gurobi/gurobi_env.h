//
// Created by Avinash on 11/30/2020.
//

#ifndef STUKA_INCLUDE_STUKA_GUROBI_GUROBI_ENV_H_
#define STUKA_INCLUDE_STUKA_GUROBI_GUROBI_ENV_H_

#include <memory>

#include <gurobi_c++.h>

namespace stuka {

/**
 * Gurobi environment.
 *
 * Creates a globally available Gurobi environment. This is meant to reduce the overhead of repeatedly creating Gurobi
 * instances.
 */
class GurobiEnvironment {
 public:
  GurobiEnvironment() {}
  ~GurobiEnvironment() {}

  GRBEnv getEnv() {
    if (!env_) env_ = std::make_unique<GRBEnv>();
    return *env_;
  }

 private:
  std::unique_ptr<GRBEnv> env_;
};

/// Globally available Gurobi environment
extern GurobiEnvironment gurobi_env;

}

#endif //STUKA_INCLUDE_STUKA_GUROBI_GUROBI_ENV_H_
