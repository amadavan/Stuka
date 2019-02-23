//
// Created by Avinash Madavan on 1/7/19.
//

#ifndef STUKA_QP_GUROBI_QP_H
#define STUKA_QP_GUROBI_QP_H

#include <algorithm>

#include <gurobi_c++.h>

#include "qp.h"
#include "base_qp.h"
#include "../constants.h"

namespace stuka { namespace QP {
  class GurobiQuadraticProgram : public BaseQuadraticProgram {
  public:
    ~GurobiQuadraticProgram() override;

    explicit GurobiQuadraticProgram(const QuadraticProgram &prog);

    void setObjective(const std::shared_ptr<Eigen::SparseMatrix<double>> &Q,
                      const std::shared_ptr<Eigen::VectorXd> &c) override;

    void setRHS(const std::shared_ptr<Eigen::VectorXd> &b_ub, const std::shared_ptr<Eigen::VectorXd> &b_eq) override;

    void setBounds(const std::shared_ptr<Eigen::VectorXd> &lb, const std::shared_ptr<Eigen::VectorXd> &ub) override;

    void addVar(double c,
        std::shared_ptr<Eigen::VectorXd> a_ub,
        std::shared_ptr<Eigen::VectorXd> a_eq,
        double lb,
        double ub) override;

    void addVars(std::shared_ptr<Eigen::VectorXd> c,
                 std::shared_ptr<Eigen::SparseMatrix<double>> A_ub,
                 std::shared_ptr<Eigen::SparseMatrix<double>> A_eq,
                 std::shared_ptr<Eigen::VectorXd> lb,
                 std::shared_ptr<Eigen::VectorXd> ub) override;

    void removeVar(size_t var) override;

    void removeVars(size_t index, size_t n_remove) override;

    void removeBackVars(size_t n_remove) override;

    void addConstr_ub(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) override;

    void addConstrs_ub(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                       const std::shared_ptr<Eigen::VectorXd> &b) override;

    void removeConstr_ub(size_t index) override;

    void removeConstrs_ub(size_t index, size_t n_remove) override;

    void addConstr_eq(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) override;

    void addConstrs_eq(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                       const std::shared_ptr<Eigen::VectorXd> &b) override;

    void removeConstr_eq(size_t index) override;

    void removeConstrs_eq(size_t index, size_t n_remove) override;

  private:
    GRBEnv env_;
    GRBModel model_;

    GRBVar *vars_;
    GRBConstr *eqconstr_;
    GRBConstr *ubconstr_;

    long n_dim_;
    long n_con_ub_;
    long n_con_eq_;

    long n_alloc_;                // Number of allocated variables (to prevent unnecessary copies)
    long n_alloc_ub_;             // Number of allocated ub constraints (to prevent unnecessary copies)
    long n_alloc_eq_;             // Number of allocated eq constraints (to prevent unnecessary copies)

    friend class GurobiSolver;
  };
}}

#endif //STUKA_QP_GUROBI_QP_H
