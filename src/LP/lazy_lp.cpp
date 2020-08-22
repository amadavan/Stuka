//
// Created by Avinash on 8/19/2020.
//

#include <stuka/lp/lazy_lp.h>

stuka::LP::LazyLinearProgram::LazyLinearProgram() {}

void stuka::LP::LazyLinearProgram::initialize(const stuka::LP::LinearProgram &prog) {
  A_ub_ = prog.A_ub;
  b_ub_ = prog.b_ub;

  n_ub_ = (b_ub_) ? b_ub_->size() : 0;
  n_active_ = 0;
  ub_activity_ = Eigen::Matrix<bool, Eigen::Dynamic, 1>::Constant(n_ub_, false);
  ub_index_ = Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Zero(n_ub_);
}

void stuka::LP::LazyLinearProgram::setObjective(const std::shared_ptr<Eigen::VectorXd> &c) {
  prog_->setObjective(c);
}

void stuka::LP::LazyLinearProgram::setRHS(const std::shared_ptr<Eigen::VectorXd> &b_ub,
                                          const std::shared_ptr<Eigen::VectorXd> &b_eq) {
  std::shared_ptr<Eigen::VectorXd> b_ub_active = nullptr;
  if (b_ub) {
    assert(b_ub_->size() == b_ub->size());
    b_ub_ = b_ub;
    b_ub_active = std::make_shared<Eigen::VectorXd>(n_active_);
    for (size_t i = 0; i < b_ub_->size(); ++i)
      if (ub_activity_.coeff(i))
        b_ub_active->coeffRef(ub_index_.coeff(i)) = b_ub_->coeff(i);
  }

  prog_->setRHS(b_ub_active, b_eq);
}

void stuka::LP::LazyLinearProgram::setBounds(const std::shared_ptr<Eigen::VectorXd> &lb,
                                             const std::shared_ptr<Eigen::VectorXd> &ub) {
  prog_->setBounds(lb, ub);
}

void stuka::LP::LazyLinearProgram::addVar(double c,
                                          std::shared_ptr<Eigen::VectorXd> a_ub,
                                          std::shared_ptr<Eigen::VectorXd> a_eq,
                                          double lb,
                                          double ub) {
  A_ub_->conservativeResize(A_ub_->rows(), A_ub_->cols() + 1);
  for (size_t i = 0; i < A_ub_->rows(); ++i)
    if (abs(a_ub->coeff(i)) > 1e-16)
      A_ub_->coeffRef(A_ub_->cols() - 1, i) = a_ub->coeff(i);

  std::shared_ptr<Eigen::VectorXd> _a_ub = nullptr;
  if (a_ub) {
    _a_ub = std::make_shared<Eigen::VectorXd>(n_active_);
    for (size_t i = 0; i < a_ub->size(); ++i)
      if (ub_activity_.coeff(i) && abs(a_ub->coeff(i)) > 1e-16)
        _a_ub->coeffRef(ub_index_.coeff(i)) = a_ub->coeff(i);
  }

  prog_->addVar(c, _a_ub, a_eq, lb, ub);
}

void stuka::LP::LazyLinearProgram::addVars(std::shared_ptr<Eigen::VectorXd> c,
                                           std::shared_ptr<Eigen::SparseMatrix<double>> A_ub,
                                           std::shared_ptr<Eigen::SparseMatrix<double>> A_eq,
                                           std::shared_ptr<Eigen::VectorXd> lb,
                                           std::shared_ptr<Eigen::VectorXd> ub) {
  size_t n_add =
      (c) ? c->size() : (lb) ? lb->size() : (ub) ? ub->size() : (A_ub) ? A_ub->cols() : (A_eq) ? A_eq->cols() : 0;

  A_ub_->conservativeResize(A_ub_->rows(), A_ub_->cols() + n_add);

  std::shared_ptr<Eigen::SparseMatrix<double>> A_ub_active = nullptr;
  if (A_ub) {
    A_ub_active = std::make_shared<Eigen::SparseMatrix<double>>(n_active_, n_add);
    for (size_t i = 0; i < n_add; ++i) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub, i); it; ++it) {
        A_ub_->coeffRef(it.row(), A_ub_->cols() - n_add + i) = it.value();
        if (ub_activity_.coeff(i))
          A_ub_active->coeffRef(ub_index_.coeff(it.row()), i) = it.value();
      }
    }
  }

  prog_->addVars(c, A_ub_active, A_eq, lb, ub);
}

void stuka::LP::LazyLinearProgram::removeVar(size_t var) {
  removeVars(var, 1);
}

void stuka::LP::LazyLinearProgram::removeVars(size_t index, size_t n_remove) {
  if (A_ub_) {
    size_t n_dim_orig = A_ub_->cols();
    std::shared_ptr<Eigen::SparseMatrix<double>>
        A_ub = std::make_shared<Eigen::SparseMatrix<double>>(A_ub_->rows(), n_dim_orig - n_remove);
    for (size_t i = 0; i < index; ++i) {
      A_ub->startVec(i);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub_, i); it; ++it)
        A_ub->insertBack(it.row(), i) = it.value();
    }
    for (size_t i = index; i < n_dim_orig - n_remove; ++i) {
      A_ub->startVec(i);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub_, i + n_remove); it; ++it)
        A_ub->insertBack(it.row(), i) = it.value();
    }
    A_ub->finalize();

    A_ub_ = A_ub;
  }

  prog_->removeVars(index, n_remove);
}

void stuka::LP::LazyLinearProgram::removeBackVars(size_t n_remove) {
  if (A_ub_) {
    size_t n_dim_orig = A_ub_->cols();
    std::shared_ptr<Eigen::SparseMatrix<double>>
        A_ub = std::make_shared<Eigen::SparseMatrix<double>>(A_ub_->rows(), n_dim_orig - n_remove);
    for (size_t i = 0; i < n_dim_orig - n_remove; ++i) {
      A_ub->startVec(i);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub_, i); it; ++it)
        A_ub->insertBack(it.row(), i) = it.value();
    }
    A_ub->finalize();

    A_ub_ = A_ub;
  }

  prog_->removeBackVars(n_remove);
}

void stuka::LP::LazyLinearProgram::addConstr_ub(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) {
  const size_t original_size = n_ub_;
  const size_t n_add = 1;

  if (A_ub_) {
    A_ub_->conservativeResize(original_size + n_add, A_ub_->cols());
    b_ub_->conservativeResize(original_size + n_add);
  }
  else {
    A_ub_ = std::make_shared<Eigen::SparseMatrix<double>>(1, a->size());
    b_ub_ = std::make_shared<Eigen::VectorXd>(1);
  }

  for (size_t i = 0; i < a->size(); ++i)
    if (abs(a->coeff(i)) > 1e-16)
      A_ub_->coeffRef(original_size, i) = a->coeff(i);
  b_ub_->coeffRef(original_size) = b;

  ub_activity_.conservativeResize(original_size + n_add);
  ub_activity_.coeffRef(original_size + 1) = false;
  ub_index_.conservativeResize(original_size + n_add);

  n_ub_ += n_add;
}

void stuka::LP::LazyLinearProgram::addConstrs_ub(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                                                 const std::shared_ptr<Eigen::VectorXd> &b) {
  const size_t original_size = n_ub_;
  const size_t n_add = b->size();

  if (A_ub_) {
    A_ub_->conservativeResize(original_size + n_add, A_ub_->cols());
    b_ub_->conservativeResize(original_size + n_add);

    for (size_t i = 0; i < A->cols(); ++i)
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A, i); it; ++it)
        A_ub_->coeffRef(original_size + it.row(), i) = it.value();
    b_ub_->tail(n_add) = *b;
  } else {
    A_ub_ = std::make_shared<Eigen::SparseMatrix<double>>(*A);
    b_ub_ = std::make_shared<Eigen::VectorXd>(*b);
  }

  ub_activity_.conservativeResize(original_size + n_add);
  for (size_t i = 0; i < n_add; ++i)
    ub_activity_.coeffRef(original_size + i) = false;
  ub_index_.conservativeResize(original_size + n_add);

  n_ub_ += n_add;
}

void stuka::LP::LazyLinearProgram::removeConstr_ub(size_t index) {
  removeConstrs_ub(index, 1);
}

void stuka::LP::LazyLinearProgram::removeConstrs_ub(size_t index, size_t n_remove) {
  std::shared_ptr<Eigen::SparseMatrix<double>>
      A_ub = std::make_shared<Eigen::SparseMatrix<double>>(A_ub_->rows() - n_remove, A_ub_->cols());
  std::shared_ptr<Eigen::VectorXd> b_ub = std::make_shared<Eigen::VectorXd>(b_ub_->size() - n_remove);

  for (size_t i = 0; i < A_ub_->cols(); ++i) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub_, i); it; ++it) {
      if (it.row() < index) A_ub->insertBack(it.row(), i) = it.value();
      else if (it.row() >= index + n_remove) A_ub->insertBack(it.row() - n_remove, i) = it.value();
    }
  }
  b_ub->head(index) = b_ub_->head(index);
  b_ub->tail(b_ub->size() - index) = b_ub_->tail(b_ub_->size() - index - n_remove);

  A_ub_ = A_ub;
  b_ub_ = b_ub;

  size_t n_deactivated = 0;

  for (size_t i = 0; i < n_remove; ++i) {
    if (ub_activity_.coeff(index + i)) {
      prog_->removeConstr_ub(ub_index_.coeff(index + i)); // Remove constraint
      n_deactivated++;

      // Decrement indices
      for (size_t j = 0; j < ub_index_.size(); ++j)
        if (ub_activity_.coeff(j) && ub_index_.coeff(j) > ub_index_.coeff(i))
          ub_index_.coeffRef(j) -= 1;
    }
  }

  Eigen::Matrix<bool, Eigen::Dynamic, 1> ub_activity(n_ub_ - n_remove);
  Eigen::Matrix<size_t, Eigen::Dynamic, 1> ub_index(n_ub_ - n_remove);
  ub_activity.head(index) = ub_activity_.head(index);
  ub_index.head(index) = ub_index_.head(index);
  ub_activity.tail(n_ub_ - index - n_remove) = ub_activity_.tail(n_ub_ - index - n_remove);
  ub_index.tail(ub_index_.size() - index - n_remove) = ub_index_.tail(ub_index_.size() - index - n_remove);
  ub_activity_ = ub_activity;
  ub_index_ = ub_index;

  n_active_ -= n_deactivated;
  n_ub_ -= n_remove;
}

void stuka::LP::LazyLinearProgram::addConstr_eq(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) {
  prog_->addConstr_eq(a, b);
}

void stuka::LP::LazyLinearProgram::addConstrs_eq(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                                                 const std::shared_ptr<Eigen::VectorXd> &b) {
  prog_->addConstrs_eq(A, b);
}

void stuka::LP::LazyLinearProgram::removeConstr_eq(size_t index) {
  prog_->removeConstr_eq(index);
}

void stuka::LP::LazyLinearProgram::removeConstrs_eq(size_t index, size_t n_remove) {
  prog_->removeConstrs_eq(index, n_remove);
}

Eigen::VectorXd stuka::LP::LazyLinearProgram::convertState(const Eigen::VectorXd &x) {
  return x;
}

Eigen::VectorXd stuka::LP::LazyLinearProgram::revertState(const Eigen::VectorXd &x) {
  return x;
}

void stuka::LP::LazyLinearProgram::setLP(stuka::LP::BaseLinearProgram *prog) {
  prog_ = prog;
}

bool stuka::LP::LazyLinearProgram::isActive(size_t index) {
  return ub_activity_.coeff(index);
}

void stuka::LP::LazyLinearProgram::setConstrActive(size_t index) {
  if (ub_activity_.coeff(index)) return;

  std::shared_ptr<Eigen::VectorXd> a = std::make_shared<Eigen::VectorXd>(A_ub_->row(index));
  prog_->addConstr_ub(a, b_ub_->coeff(index));
  ub_activity_.coeffRef(index) = true;
  ub_index_.coeffRef(index) = n_active_++;
}

void stuka::LP::LazyLinearProgram::setConstrsActive(size_t index, size_t n_active) {

  throw std::runtime_error("LP::LazyLinearProgram::setConstrsActive not implemented.");
  // TODO: Add check for whether constraints are already active.
  for (size_t i = 0; i < n_active; ++i)
    if (ub_activity_.coeff(index + i))
      throw std::runtime_error(
          "LP::LazyLinearProgram::setConstrsActive not implemented with existing active constraints.");

  std::shared_ptr<Eigen::SparseMatrix<double>>
      A = std::make_shared<Eigen::SparseMatrix<double>>(A_ub_->middleRows(index, n_active));
  std::shared_ptr<Eigen::VectorXd> b = std::make_shared<Eigen::VectorXd>(b_ub_->segment(index, n_active));
  prog_->addConstrs_ub(A, b);
  for (size_t i = 0; i < n_active; ++i) {
    ub_activity_.coeffRef(index + i) = true;
    ub_index_.coeffRef(index + i) = n_active_++;
  }
}

void stuka::LP::LazyLinearProgram::setConstrsActive(Eigen::Matrix<bool, Eigen::Dynamic, 1> activity) {
  assert(n_ub_ == activity.size());

  for (size_t i = 0; i < activity.size(); ++i)
    if (activity.coeff(i) && !ub_activity_.coeff(i)) setConstrActive(i);
    else if (!activity.coeff(i) && ub_activity_.coeff(i))
      std::cout << "WARNING: setting active constraints inactive is not supported." << std::endl;
    else {}

  ub_activity_ = activity;
}

void stuka::LP::LazyLinearProgram::setConstrsActive(std::vector<size_t> activity) {
  for (size_t i = 0; i < activity.size(); ++i)
    setConstrActive(activity[i]);
}

bool stuka::LP::LazyLinearProgram::isViolated(const Eigen::VectorXd &x) {
  return x.size() == 0 || ((*A_ub_ * x - *b_ub_).array() > 0).any();
}

bool stuka::LP::LazyLinearProgram::addViolatedConstraints(const Eigen::VectorXd &x) {
  if (!A_ub_ || x.size() == 0) return false;
  Eigen::Matrix<bool, Eigen::Dynamic, 1> violated = ((*A_ub_ * x - *b_ub_).array() > 0);
  for (size_t i = 0; i < violated.size(); ++i) {
    if (violated.coeff(i) && !ub_activity_.coeff(i))
      setConstrActive(i);
  }
  return violated.any();
}

bool stuka::LP::LazyLinearProgram::addRandomConstraint() {
  for (size_t i = 0; i < n_ub_; ++i) {
    if (!ub_activity_.coeff(i)) {
      setConstrActive(i);
      return true;
    }
  }
  return false;
}

Eigen::VectorXd stuka::LP::LazyLinearProgram::getDualUB(const Eigen::VectorXd &active_dual) {
  Eigen::VectorXd dual = Eigen::VectorXd::Constant(b_ub_->size(), 0);
  for (size_t i = 0; i < n_ub_; ++i)
    if (ub_activity_.coeff(i))
      dual.coeffRef(i) = active_dual.coeff(ub_index_.coeff(i));
  return dual;
}
