//
// Created by Avinash Madavan on 2019-01-29.
//

#ifndef STUKA_PYTHON_DLP_H
#define STUKA_PYTHON_DLP_H

namespace stuka { namespace dLP {

  class PyDecomposedLinearProgram : public DecomposedLinearProgram {
  public:
    using DecomposedLinearProgram::DecomposedLinearProgram;

    PyDecomposedLinearProgram(std::vector <Eigen::VectorXd> _c,
                              std::vector <Eigen::SparseMatrix<double>> _A_ub,
                              std::vector <Eigen::VectorXd> _b_ub,
                              std::vector <Eigen::SparseMatrix<double>> _C_ub,
                              std::vector <Eigen::SparseMatrix<double>> _A_eq,
                              std::vector <Eigen::VectorXd> _b_eq,
                              std::vector <Eigen::SparseMatrix<double>> _C_eq,
                              std::vector <Eigen::VectorXd> _lb,
                              std::vector <Eigen::VectorXd> _ub) {
      size_t n = _c.size();
      size_t n_sub = n - 1;

      c = std::vector < std::shared_ptr < Eigen::VectorXd >> (n);
      A_ub = std::vector < std::shared_ptr < Eigen::SparseMatrix < double >> > (n);
      b_ub = std::vector < std::shared_ptr < Eigen::VectorXd >> (n);
      C_ub = std::vector < std::shared_ptr < Eigen::SparseMatrix < double >> > (n_sub);
      A_eq = std::vector < std::shared_ptr < Eigen::SparseMatrix < double >> > (n);
      b_eq = std::vector < std::shared_ptr < Eigen::VectorXd >> (n);
      C_eq = std::vector < std::shared_ptr < Eigen::SparseMatrix < double >> > (n_sub);
      lb = std::vector < std::shared_ptr < Eigen::VectorXd >> (n);
      ub = std::vector < std::shared_ptr < Eigen::VectorXd >> (n);

      for (size_t i = 0; i < n_sub; ++i) {
        c[i] = std::make_shared<Eigen::VectorXd>(_c[i]);
        A_ub[i] = std::make_shared < Eigen::SparseMatrix < double >> (_A_ub[i]);
        b_ub[i] = std::make_shared<Eigen::VectorXd>(_b_ub[i]);
        C_ub[i] = std::make_shared < Eigen::SparseMatrix < double >> (_C_ub[i]);
        A_eq[i] = std::make_shared < Eigen::SparseMatrix < double >> (_A_eq[i]);
        b_eq[i] = std::make_shared<Eigen::VectorXd>(_b_eq[i]);
        C_eq[i] = std::make_shared < Eigen::SparseMatrix < double >> (_C_eq[i]);
        lb[i] = std::make_shared<Eigen::VectorXd>(_lb[i]);
        ub[i] = std::make_shared<Eigen::VectorXd>(_ub[i]);
      }

      c[n_sub] = std::make_shared<Eigen::VectorXd>(_c[n_sub]);
      A_ub[n_sub] = std::make_shared < Eigen::SparseMatrix < double >> (_A_ub[n_sub]);
      b_ub[n_sub] = std::make_shared<Eigen::VectorXd>(_b_ub[n_sub]);
      A_eq[n_sub] = std::make_shared < Eigen::SparseMatrix < double >> (_A_eq[n_sub]);
      b_eq[n_sub] = std::make_shared<Eigen::VectorXd>(_b_eq[n_sub]);
      lb[n_sub] = std::make_shared<Eigen::VectorXd>(_lb[n_sub]);
      ub[n_sub] = std::make_shared<Eigen::VectorXd>(_ub[n_sub]);
    }

  };

}}

#endif //STUKA_PYTHON_DLP_H
