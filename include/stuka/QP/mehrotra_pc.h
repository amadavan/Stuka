#ifndef STUKA_QP_MEHROTRA_PC_CPP_
#define STUKA_QP_MEHROTRA_PC_CPP_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/CholmodSupport>
#include <stuka/QP/base_solver.h>
#include <stuka/LP/base_solver.h>

#include <stuka/LP/lp.h>

namespace stuka::QP {

  class MehrotraPC : public LP::BaseLPSolver {
  public:
    MehrotraPC(const LP::LinearProgram &lp, const Options &opts);
    ~MehrotraPC() = default;
  
    void iterate() override;

    bool terminate() override;

    const OptimizeState getState() override;

    LP::BaseLinearProgram &getLP() override;
  private:    
    /**
     * @brief The IPM step.
     * 
     * Specifies the step to be taken for each of the variables.
     * 
     */
    struct Step {
      Eigen::VectorXd dx;
      Eigen::VectorXd dy;
      Eigen::VectorXd dz;
      Eigen::VectorXd ds;
      Eigen::VectorXd dw;
    };

    /**
     * @brief Step size to take for primal/dual variables.
     * 
     * A simple struct describing the step size to take for primal variables
     * (x, z) and dual variables (y, s, w).
     */
    struct StepSize {
      double primal;
      double dual;
    };
    
    /**
     * @brief Right-hand-side of the linear system.
     * 
     * Specifies each of the components of the right-hand-side of the linear
     * system describing the KKT conditions. For example,
     *    Rc = c - A.T y - s + w,
     *    Rb = b - A x,
     *    Ru = u - x - z,
     *    Rxs = -XSe,
     *    Rzw = -ZWe
     */
    struct RHS {
      Eigen::VectorXd Rc;         ///< Error for gradient condition
      Eigen::VectorXd Rb;         ///< Error for feasibility condition (Ax = b)
      Eigen::VectorXd Ru;         ///< Error for upper-bound slack equality
      Eigen::VectorXd Rxs;        ///< Error for complementary slack on x/s
      Eigen::VectorXd Rzw;        ///< Error for complementary slack on z/w
    };

    /**
     * @brief Set the state to an initial interior point.
     * 
     */
    void computeInitialPoint();
    
    /**
     * @brief Solve the IPM KKT system of equations.
     * 
     * Computes a single iteration of the IPM step for a given right-hand-side.
     * This solves the system of equations specified by the KKT conditions with
     * a fixed error vector.
     * 
     * The KKT system that is solved should be the gradient condition, whose
     * error is denoted Rc, the feasibility condition, with error denoted by Rb,
     * the upper-bound slack equality constraint, with error denoted by Ru, the
     * complementary slackness conditions (with log-barrier), with errors 
     * denoted by Rxs and Rzw for the standard variables and compact upper-bound
     * variables respectively. 
     * 
     * @param rhs A vector of errors for each of the KKT condition equations.
     * @param resolve Whether the KKT system should be resolved.
     * @return Step The step to be taken for each of the variables.
     */
    Step solveKKT(RHS rhs, bool resolve = true);

    /**
     * @brief Determine an appropriate step-size for the current iterate.
     * 
     * Performs a line search to evaluate the maximal step-size, less than 1,
     * that ensures that the point remains feasible. This is step can be scaled
     * by the parameter to ensure strict feasibility.
     * 
     * @param step The step for which to evaluate the step size.
     * @param scale Scaling parameter for step size (should be <= 1).
     * @return StepSize Appropriate step-size for given step.
     */
    StepSize computeStepSize(Step step, double scale = 1.);

    /**
     * @brief Type of variable bounds.
     * 
     * A descriptor for the variable depending on the description of the bounds
     * of the variable.
     * 
     */
    enum VariableType {
      STANDARD,   ///< Lower-bound of 0 and no upper-bound
      FREE,       ///< No bounds
      COMPACT,    ///< Lower and upper bounds
      LOWER,      ///< Lower bounded (non-zero) and no upper bound
      UPPER       ///< Upper bounded and no lower-bound
    };

    /// Map of variable index (in x) to variable type
    std::map<int, VariableType> variable_types_;

    /* 
     * Problem description in the form
     *      min  c.T x,
     *      s.t. A x == b,
     *           0 <= x <= u.
     */

    /// Matrix of constraint coefficients
    Eigen::SparseMatrix<double, Eigen::RowMajor> A_;
    Eigen::VectorXd b_;             ///< Vector of constraint bounds
    Eigen::VectorXd c_;             ///< Vector of cost coefficients
    Eigen::SparseMatrix<double> Q_; ///< Matrix of quadratic costs
    Eigen::VectorXd u_;             ///< Vector of apparent bounds of compact vars
    Eigen::VectorXd bound_;         ///< Vector of actual bound shifts

    /// Map of compact variable index (in z) to index (in x)
    std::map<int, int> z_index_;

    int n_;                       ///< Number of variables
    int n_base_;                  ///< Number of variables in original problem
    int n_slack_;                 ///< Number of slack variables
    int n_compact_;               ///< Number of compact variables
    int n_free_;                  ///< Number of free variables
    int m_eq_;                    ///< Number of equality constraints
    int m_ub_;                    ///< Number of inequality cons (upper-bounded)

    double c_norm_;               ///< Norm of cost coefficients
    double b_norm_;               ///< Norm of constraint bounds
    double u_norm_;               ///< Norm of apparent bounds

    /*
     * Solver parameters 
     */
    Eigen::VectorXd x_;           ///< Primal variable
    Eigen::VectorXd y_;           ///< Dual variable
    Eigen::VectorXd z_;           ///< Upper-bound slack for compact vars
    Eigen::VectorXd s_;           ///< Lower-bound dual multiplier for x
    Eigen::VectorXd w_;           ///< Lower-bound dual multiplier for z

    double mu_;                   ///< Barrier parameter
    RHS rhs_;                     ///< Current right-hand-side vector

    /// Linear system solver for (scaled) A * At system
    Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double, Eigen::RowMajor>> ADA_t;

    double error_;                ///< Current solver error

    const Options &opts_;          ///< Solver options
    const LP::LinearProgram &prog_;   // TODO: remove this and fix the framework structure
  };
  
}

#endif // STUKA_QP_MEHROTRA_PC_CPP_