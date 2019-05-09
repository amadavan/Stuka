//
// Created by Avinash Madavan on 1/10/19.
//

#ifndef STUKA_BASE_SOLVER_H
#define STUKA_BASE_SOLVER_H

#ifndef DEFAULT_MAX_ITER
#define DEFAULT_MAX_ITER 100
#endif

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "options.h"
#include "optimize_state.h"


namespace stuka {

  class BaseSolver {
  public:
    explicit BaseSolver(const Options &opts) : opts_(opts), n_max_iter_(opts.max_iter) {
      if (opts.max_iter == 0)
        n_max_iter_ = DEFAULT_MAX_ITER;
    }

    virtual ~BaseSolver() = default;

    /* Perform one solver iteration.
     *
     * The scheme defining how to update the solver state given the current state.
     */
    virtual void iterate() = 0;

    /* Evaluate termination criterion.
     *
     * Check for convergence or divergence of the algorithm given the current state of the optimizer.
     */
    virtual bool terminate() = 0;

    /* Get the current optimizer state.
     *
     * Returns an OptimizeState structure containing solver-provided information about the current iterate.
     */
    virtual const OptimizeState getState() = 0;

    /* Solve the optimization problem.
     *
     * Iterate until the termination criterion is met or the maximum number of iterations has been reached.
     */
    virtual const OptimizeState solve() {
      unsigned int nit = 0;
      do {
        iterate();
        if (opts_.callback) {
          OptimizeState state = this->getState();
          state.nit = nit;
          opts_.callback->callback(state);
        }
      } while (++nit < n_max_iter_ && !terminate());

      OptimizeState res = getState();
      res.nit = nit;

      return res;
    }

  private:
    size_t n_max_iter_;
    const Options &opts_;
  };

}

#endif //STUKA_BASE_SOLVER_H
