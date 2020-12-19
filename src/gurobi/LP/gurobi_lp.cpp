//
// Created by Avinash Madavan on 1/3/19.
//

#include <stuka/gurobi/LP/gurobi_lp.h>

stuka::LP::GurobiLinearProgram::~GurobiLinearProgram() {
  delete[] vars_;
  if (n_con_eq_ > 0) delete[] eqconstr_;
  if (n_con_ub_ > 0) delete[] ubconstr_;
}

// TODO: use asserts to check for invalid inputs
stuka::LP::GurobiLinearProgram::GurobiLinearProgram() : model_(gurobi_env.getEnv()) {}

void stuka::LP::GurobiLinearProgram::initialize(const stuka::LP::LinearProgram &prog) {
  n_dim_ = prog.c->size();
  n_con_ub_ = prog.b_ub ? prog.b_ub->size() : 0;
  n_con_eq_ = prog.b_eq ? prog.b_eq->size() : 0;

  n_alloc_ = n_dim_;
  n_alloc_ub_ = n_con_ub_;
  n_alloc_eq_ = n_con_eq_;

  // Set bounds
  Eigen::VectorXd lb = Eigen::VectorXd::Constant(n_dim_, -GRB_INFINITY);
  Eigen::VectorXd ub = Eigen::VectorXd::Constant(n_dim_, GRB_INFINITY);

  if (prog.lb)
    for (size_t i = 0; i < n_dim_; ++i)
      lb.coeffRef(i) = (prog.lb->coeff(i) == -INF) ? -GRB_INFINITY : prog.lb->coeff(i);

  if (prog.ub)
    for (size_t i = 0; i < n_dim_; ++i)
      ub.coeffRef(i) = (prog.ub->coeff(i) == INF) ? GRB_INFINITY : prog.ub->coeff(i);

  // Create variables
  vars_ = model_.addVars(lb.data(), ub.data(), prog.c->data(), nullptr, nullptr, n_dim_);

  // Create inequality constraints
  if (n_con_ub_ > 0) {
    ubconstr_ = model_.addConstrs((int) n_con_ub_);
    for (size_t i = 0; i < n_con_ub_; ++i) {
      (ubconstr_ + i)->set(GRB_CharAttr_Sense, '<');
      (ubconstr_ + i)->set(GRB_DoubleAttr_RHS, prog.b_ub->coeff(i));
    }
    for (size_t i = 0; i < n_dim_; ++i) {
      for (Eigen::SparseMatrix<double>::InnerIterator it_ub(*prog.A_ub, i); it_ub; ++it_ub)
        model_.chgCoeff(*(ubconstr_ + it_ub.row()), *(vars_ + i), it_ub.value());
    }
  } else {
    ubconstr_ = nullptr;
  }

  // Create equality constraints
  if (n_con_eq_ > 0) {
    eqconstr_ = model_.addConstrs((int) n_con_eq_);
    for (size_t i = 0; i < n_con_eq_; ++i) {
      (eqconstr_ + i)->set(GRB_CharAttr_Sense, '=');
      (eqconstr_ + i)->set(GRB_DoubleAttr_RHS, prog.b_eq->coeff(i));
    }
    for (size_t i = 0; i < n_dim_; ++i) {
      for (Eigen::SparseMatrix<double>::InnerIterator it_eq(*prog.A_eq, i); it_eq; ++it_eq)
        model_.chgCoeff(*(eqconstr_ + it_eq.row()), *(vars_ + i), it_eq.value());
    }
  } else {
    eqconstr_ = nullptr;
  }

  model_.set(GRB_IntParam_OutputFlag, false);
  model_.set(GRB_DoubleParam_OptimalityTol, GUROBI_TOLERANCE);
  model_.set(GRB_DoubleParam_FeasibilityTol, GUROBI_TOLERANCE);
}

void stuka::LP::GurobiLinearProgram::setObjective(const std::unique_ptr<Eigen::VectorXd> &c) {
  // Set objective
  GRBLinExpr obj;
  obj.addTerms(c->data(), vars_, (int) n_dim_);
  model_.setObjective(obj);

}

void stuka::LP::GurobiLinearProgram::setRHS(const std::unique_ptr<Eigen::VectorXd> &b_ub,
                                            const std::unique_ptr<Eigen::VectorXd> &b_eq) {
  if (b_ub && n_con_ub_ > 0)
    for (size_t i = 0; i < n_con_ub_; ++i)
      (ubconstr_ + i)->set(GRB_DoubleAttr_RHS, b_ub->coeff(i));

  if (b_eq && n_con_eq_ > 0)
    for (size_t i = 0; i < n_con_eq_; ++i)
      (eqconstr_ + i)->set(GRB_DoubleAttr_RHS, b_eq->coeff(i));
}

void stuka::LP::GurobiLinearProgram::setBounds(const std::unique_ptr<Eigen::VectorXd> &lb,
                                               const std::unique_ptr<Eigen::VectorXd> &ub) {
  if (lb) model_.set(GRB_DoubleAttr_LB, vars_, lb->data(), (int) n_dim_);
  if (ub) model_.set(GRB_DoubleAttr_LB, vars_, ub->data(), (int) n_dim_);
}

void stuka::LP::GurobiLinearProgram::addVar(double c,
                                            const std::unique_ptr<Eigen::VectorXd> &a_ub,
                                            const std::unique_ptr<Eigen::VectorXd> &a_eq,
                                            double lb,
                                            double ub) {

  if (n_dim_ + 1 > n_alloc_) {
    n_alloc_ = n_dim_ + 1;
    GRBVar *vars = new GRBVar[n_alloc_];
    std::copy(vars_, vars_ + n_dim_, vars);
    delete[] vars_;
    vars_ = vars;
  }

  vars_[n_dim_] = model_.addVar((lb == -INF) ? -GRB_INFINITY : lb,
                                (ub == INF) ? GRB_INFINITY : ub,
                                c, GRB_CONTINUOUS);

  if (n_con_ub_ > 0 && a_ub) {
    for (size_t i = 0; i < n_con_ub_; ++i)
      if (abs(a_ub->coeff(i)) > 1e-8)
        model_.chgCoeff(ubconstr_[i], vars_[n_dim_], a_ub->coeff(i));
  }

  if (n_con_eq_ > 0 && a_eq) {
    for (size_t i = 0; i < n_con_eq_; ++i)
      if (abs(a_eq->coeff(i)) > 1e-8)
        model_.chgCoeff(eqconstr_[i], vars_[n_dim_], a_eq->coeff(i));
  }

  n_dim_ += 1;

}

void stuka::LP::GurobiLinearProgram::addVars(const std::unique_ptr<Eigen::VectorXd> &c,
                                             const std::unique_ptr<Eigen::SparseMatrix<double>> &A_ub,
                                             const std::unique_ptr<Eigen::SparseMatrix<double>> &A_eq,
                                             const std::unique_ptr<Eigen::VectorXd> &lb_,
                                             const std::unique_ptr<Eigen::VectorXd> &ub_) {

  size_t n_add = (c) ? c->size() : (lb_) ? lb_->size() : (ub_) ? ub_->size() :
                                                         (A_ub) ? A_ub->cols() : (A_eq) ? A_eq->cols() : 0;
  if (n_add == 0) return;

  if (n_dim_ + n_add > n_alloc_) {
    n_alloc_ = n_dim_ + n_add;
    GRBVar *vars = new GRBVar[n_alloc_];
    std::copy(vars_, vars_ + n_dim_, vars);
    delete[] vars_;
    vars_ = vars;
  }

  // Set bounds
  Eigen::VectorXd lb = Eigen::VectorXd::Constant(n_add, -GRB_INFINITY);
  Eigen::VectorXd ub = Eigen::VectorXd::Constant(n_add, GRB_INFINITY);

  if (lb_)
    for (size_t i = 0; i < n_add; ++i)
      lb.coeffRef(i) = (lb_->coeff(i) == -INF) ? -GRB_INFINITY : lb_->coeff(i);

  if (ub_)
    for (size_t i = 0; i < n_add; ++i)
      ub.coeffRef(i) = (ub_->coeff(i) == INF) ? GRB_INFINITY : ub_->coeff(i);

  // Add variables
  GRBVar *vars_add = model_.addVars(lb.data(), ub.data(), (c) ? c->data() : nullptr, nullptr, nullptr, n_add);
  std::copy(vars_add, vars_add + n_add, vars_ + n_dim_);
  delete[] vars_add;

  // Add constraints
  if (n_con_ub_ > 0 && A_ub) {
    for (size_t i = 0; i < n_add; ++i) {
      for (Eigen::SparseMatrix<double>::InnerIterator it_ub(*A_ub, i); it_ub; ++it_ub)
        model_.chgCoeff(*(ubconstr_ + it_ub.row()), *(vars_ + n_dim_ + i), it_ub.value());
    }
  }

  if (n_con_eq_ > 0 && A_eq) {
    for (size_t i = 0; i < n_add; ++i) {
      for (Eigen::SparseMatrix<double>::InnerIterator it_eq(*A_eq, i); it_eq; ++it_eq)
        model_.chgCoeff(*(eqconstr_ + it_eq.row()), *(vars_ + n_dim_ + i), it_eq.value());
    }
  }

  n_dim_ += n_add;

}

void stuka::LP::GurobiLinearProgram::removeVar(const size_t index) {
  model_.remove(vars_[index]);
  std::copy(vars_ + index + 1, vars_ + n_dim_, vars_ + index);
  n_dim_ -= 1;
}

void stuka::LP::GurobiLinearProgram::removeVars(const size_t index, const size_t n_remove) {

  for (size_t i = 0; i < n_remove; ++i)
    model_.remove(vars_[index + i]);
  std::copy(vars_ + index + n_remove, vars_ + n_dim_, vars_ + index);

  n_dim_ -= n_remove;

}

void stuka::LP::GurobiLinearProgram::removeBackVars(const size_t n_remove) {
  for (size_t i = 0; i < n_remove; ++i)
    model_.remove(vars_[n_dim_ - i - 1]);

  n_dim_ -= n_remove;
}

void stuka::LP::GurobiLinearProgram::addConstr_ub(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) {

  if (n_con_ub_ + 1 > n_alloc_ub_) {
    n_alloc_ub_ = n_con_ub_ + 1;
    GRBConstr *cons = new GRBConstr[n_alloc_ub_];
    std::copy(ubconstr_, ubconstr_ + n_con_ub_, cons);
    delete[] ubconstr_;
    ubconstr_ = cons;
  }

  GRBLinExpr expr;
  for (size_t i = 0; i < n_dim_; ++i)
    if (abs(a->coeff(i)) > 1e-8)
      expr += a->coeff(i) * vars_[i];
  ubconstr_[n_con_ub_] = model_.addConstr(expr, '<', b);

  n_con_ub_ += 1;

}

void stuka::LP::GurobiLinearProgram::addConstrs_ub(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                                                   const std::unique_ptr<Eigen::VectorXd> &b) {

  size_t n_add = b->size();

  if (n_con_ub_ + n_add > n_alloc_ub_) {
    n_alloc_ub_ = n_con_ub_ + n_add;
    GRBConstr *cons = new GRBConstr[n_alloc_ub_];
    std::copy(ubconstr_, ubconstr_ + n_con_ub_, cons);
    delete[] ubconstr_;
    ubconstr_ = cons;
  }

  GRBConstr *newcons = model_.addConstrs(n_add);
  std::copy(newcons, newcons + n_add, ubconstr_ + n_con_ub_);
  delete[] newcons;

  for (size_t i = 0; i < n_add; ++i) {
    (ubconstr_ + n_con_ub_ + i)->set(GRB_CharAttr_Sense, '<');
    (ubconstr_ + n_con_ub_ + i)->set(GRB_DoubleAttr_RHS, b->coeff(i));
  }
  for (size_t i = 0; i < n_dim_; ++i)
    for (Eigen::SparseMatrix<double>::InnerIterator it_ub(*A, i); it_ub; ++it_ub)
      model_.chgCoeff(*(ubconstr_ + n_con_ub_ + it_ub.row()), *(vars_ + i), it_ub.value());

  n_con_ub_ += n_add;

}
void stuka::LP::GurobiLinearProgram::removeConstr_ub(const size_t index) {
  removeConstrs_ub(index, 1);
}

void stuka::LP::GurobiLinearProgram::removeConstrs_ub(const size_t index, const size_t n_remove) {

  for (size_t i = 0; i < n_remove; ++i)
    model_.remove(ubconstr_[index + i]);
  std::copy(ubconstr_ + index + n_remove, ubconstr_ + n_con_ub_, ubconstr_ + index);

  n_con_ub_ -= n_remove;

}

void stuka::LP::GurobiLinearProgram::addConstr_eq(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) {

  if (n_con_eq_ + 1 > n_alloc_eq_) {
    n_alloc_eq_ = n_con_eq_ + 1;
    GRBConstr *cons = new GRBConstr[n_alloc_eq_];
    std::copy(eqconstr_, eqconstr_ + n_con_eq_, cons);
    delete[] eqconstr_;
    eqconstr_ = cons;
  }

  GRBLinExpr expr;
  for (size_t i = 0; i < n_dim_; ++i)
    if (abs(a->coeff(i)) > 1e-8)
      expr += a->coeff(i) * vars_[i];
  eqconstr_[n_con_eq_] = model_.addConstr(expr, '=', b);

  n_con_eq_ += 1;

}

void stuka::LP::GurobiLinearProgram::addConstrs_eq(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                                                   const std::unique_ptr<Eigen::VectorXd> &b) {

  size_t n_add = b->size();

  if (n_con_eq_ + n_add > n_alloc_eq_) {
    n_alloc_eq_ = n_con_eq_ + n_add;
    GRBConstr *cons = new GRBConstr[n_alloc_eq_];
    std::copy(eqconstr_, eqconstr_ + n_con_eq_, cons);
    delete[] eqconstr_;
    eqconstr_ = cons;
  }

  GRBConstr *newcons = model_.addConstrs(n_add);
  std::copy(newcons, newcons + n_add, eqconstr_ + n_con_eq_);
  delete[] newcons;

  for (size_t i = 0; i < n_add; ++i) {
    (eqconstr_ + n_con_eq_ + i)->set(GRB_CharAttr_Sense, '=');
    (eqconstr_ + n_con_eq_ + i)->set(GRB_DoubleAttr_RHS, b->coeff(i));
  }
  for (size_t i = 0; i < n_dim_; ++i)
    for (Eigen::SparseMatrix<double>::InnerIterator it_eq(*A, i); it_eq; ++it_eq)
      model_.chgCoeff(*(eqconstr_ + n_con_eq_ + it_eq.row()), *(vars_ + i), it_eq.value());

  n_con_eq_ += n_add;

}

void stuka::LP::GurobiLinearProgram::removeConstr_eq(const size_t index) {
  removeConstrs_eq(index, 1);
}

void stuka::LP::GurobiLinearProgram::removeConstrs_eq(const size_t index, const size_t n_remove) {

  for (size_t i = 0; i < n_remove; ++i)
    model_.remove(eqconstr_[index + i]);
  std::copy(eqconstr_ + index + n_remove, eqconstr_ + n_con_eq_, eqconstr_ + index);

  n_con_eq_ -= n_remove;

}

Eigen::VectorXd stuka::LP::GurobiLinearProgram::convertState(const Eigen::VectorXd &x) {
  return x;
}

Eigen::VectorXd stuka::LP::GurobiLinearProgram::revertState(const Eigen::VectorXd &x) {
  return x;
}
