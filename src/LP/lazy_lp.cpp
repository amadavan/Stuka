//
// Created by Avinash on 8/19/2020.
//

#include <stuka/LP/lazy_lp.h>

stuka::LP::LazyLinearProgram::LazyLinearProgram() {}

void stuka::LP::LazyLinearProgram::initialize(const stuka::LP::LinearProgram &prog) {
  A_ub_ = util::SparseOps::unique_copy(prog.A_ub);
  b_ub_ = util::DenseOps::unique_copy(prog.b_ub);

  n_dim_ = (prog.c) ? prog.c->size() : (prog.A_ub) ? prog.A_ub->cols() : (prog.A_eq) ? prog.A_eq->cols() : (prog.lb) ? prog.lb->size() : (prog.ub) ? prog.ub->size() : 0;
  n_ub_ = (b_ub_) ? b_ub_->size() : 0;
  n_active_ = 0;
  ub_activity_ = Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(n_ub_, false);
  ub_index_ = Eigen::Array<size_t, Eigen::Dynamic, 1>::Zero(n_ub_);
}

void stuka::LP::LazyLinearProgram::setObjective(const std::unique_ptr<Eigen::VectorXd> &c) {
  prog_->setObjective(c);
}

void stuka::LP::LazyLinearProgram::setRHS(const std::unique_ptr<Eigen::VectorXd> &b_ub,
                                          const std::unique_ptr<Eigen::VectorXd> &b_eq) {
  std::unique_ptr<Eigen::VectorXd> b_ub_active = nullptr;
  if (b_ub) {
    assert(b_ub_->size() == b_ub->size());
    b_ub_ = std::make_unique<Eigen::VectorXd>(*b_ub);
    b_ub_active = std::make_unique<Eigen::VectorXd>(n_active_);
    for (size_t i = 0; i < b_ub_->size(); ++i)
      if (ub_activity_.coeff(i))
        b_ub_active->coeffRef(ub_index_.coeff(i)) = b_ub_->coeff(i);
  }

  prog_->setRHS(b_ub_active, b_eq);
}

void stuka::LP::LazyLinearProgram::setBounds(const std::unique_ptr<Eigen::VectorXd> &lb,
                                             const std::unique_ptr<Eigen::VectorXd> &ub) {
  prog_->setBounds(lb, ub);
}

void stuka::LP::LazyLinearProgram::addVar(double c,
                                          const std::unique_ptr<Eigen::VectorXd> &a_ub,
                                          const std::unique_ptr<Eigen::VectorXd> &a_eq,
                                          double lb,
                                          double ub) {
  A_ub_->conservativeResize(A_ub_->rows(), n_dim_ + 1);
  for (size_t i = 0; i < A_ub_->rows(); ++i)
    if (abs(a_ub->coeff(i)) > 1e-16)
      A_ub_->coeffRef(n_dim_ - 1, i) = a_ub->coeff(i);

  std::unique_ptr<Eigen::VectorXd> _a_ub = nullptr;
  if (a_ub) {
    _a_ub = std::make_unique<Eigen::VectorXd>(n_active_);
    for (size_t i = 0; i < a_ub->size(); ++i)
      if (ub_activity_.coeff(i) && abs(a_ub->coeff(i)) > 1e-16)
        _a_ub->coeffRef(ub_index_.coeff(i)) = a_ub->coeff(i);
  }

  n_dim_++;

  prog_->addVar(c, _a_ub, a_eq, lb, ub);
}

void stuka::LP::LazyLinearProgram::addVars(const std::unique_ptr<Eigen::VectorXd> &c,
                                           const std::unique_ptr<Eigen::SparseMatrix<double>> &A_ub,
                                           const std::unique_ptr<Eigen::SparseMatrix<double>> &A_eq,
                                           const std::unique_ptr<Eigen::VectorXd> &lb,
                                           const std::unique_ptr<Eigen::VectorXd> &ub) {
  size_t n_add =
      (c) ? c->size() : (lb) ? lb->size() : (ub) ? ub->size() : (A_ub) ? A_ub->cols() : (A_eq) ? A_eq->cols() : 0;

  // A_ub_->conservativeResize(A_ub_->rows(), A_ub_->cols() + n_add);

  std::unique_ptr<Eigen::SparseMatrix<double>> _A_ub = std::make_unique<Eigen::SparseMatrix<double>>(n_ub_, n_dim_ + n_add);
  if (A_ub_) {
    for (size_t i = 0; i < n_dim_; ++i) {
      _A_ub->startVec(i);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub_, i); it; ++it)
        _A_ub->insertBack(it.row(), i) = it.value();
    }
  }

  std::unique_ptr<Eigen::SparseMatrix<double>> A_ub_active = nullptr;
  if (A_ub) {
    A_ub_active = std::make_unique<Eigen::SparseMatrix<double>>(n_active_, n_add);
    for (size_t i = 0; i < n_add; ++i) {
      _A_ub->startVec(n_dim_ + i);

      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub, i); it; ++it) {
        _A_ub->insertBack(it.row(), n_dim_ + i) = it.value();
        if (ub_activity_.coeff(it.row()))
          A_ub_active->coeffRef(ub_index_.coeff(it.row()), i) = it.value();
      }
    }
  }

  A_ub_ = std::move(_A_ub);
  n_dim_ += n_add;

  prog_->addVars(c, A_ub_active, A_eq, lb, ub);
}

void stuka::LP::LazyLinearProgram::removeVar(size_t var) {
  removeVars(var, 1);
}

void stuka::LP::LazyLinearProgram::removeVars(size_t index, size_t n_remove) {
  if (A_ub_) {
    std::unique_ptr<Eigen::SparseMatrix<double>>
        A_ub = std::make_unique<Eigen::SparseMatrix<double>>(A_ub_->rows(), n_dim_ - n_remove);
    A_ub->reserve(A_ub_->nonZeros());
    for (size_t i = 0; i < index; ++i) {
      A_ub->startVec(i);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub_, i); it; ++it)
        A_ub->insertBack(it.row(), i) = it.value();
    }
    for (size_t i = index; i < n_dim_ - n_remove; ++i) {
      A_ub->startVec(i);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub_, i + n_remove); it; ++it)
        A_ub->insertBack(it.row(), i) = it.value();
    }
    A_ub->finalize();

    A_ub_ = std::move(A_ub);
  }

  n_dim_ -= n_remove;
  prog_->removeVars(index, n_remove);
}

void stuka::LP::LazyLinearProgram::removeBackVars(size_t n_remove) {
  if (A_ub_) removeVars(n_dim_ - n_remove, n_remove);

  prog_->removeBackVars(n_remove);
}

void stuka::LP::LazyLinearProgram::addConstr_ub(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) {
  const size_t n_add = 1;

  if (A_ub_) {
    A_ub_->conservativeResize(n_ub_ + n_add, n_dim_);
    b_ub_->conservativeResize(n_ub_ + n_add);
  }
  else {
    A_ub_ = std::make_unique<Eigen::SparseMatrix<double>>(1, a->size());
    b_ub_ = std::make_unique<Eigen::VectorXd>(1);
  }

  for (size_t i = 0; i < a->size(); ++i)
    if (abs(a->coeff(i)) > 1e-14)
      A_ub_->coeffRef(n_ub_, i) = a->coeff(i);
  b_ub_->coeffRef(n_ub_) = b;

  ub_activity_.conservativeResize(n_ub_ + n_add);
  ub_activity_.coeffRef(n_ub_ + 1) = false;
  ub_index_.conservativeResize(n_ub_ + n_add);

  n_ub_ += n_add;
}

void stuka::LP::LazyLinearProgram::addConstrs_ub(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                                                 const std::unique_ptr<Eigen::VectorXd> &b) {
  const size_t n_add = b->size();

  if (A_ub_) {
    A_ub_ = std::make_unique<Eigen::SparseMatrix<double>>(util::SparseOps::vstack({*A_ub_, *A}));

//    A_ub_->conservativeResize(n_ub_ + n_add, A_ub_->cols());
//    for (size_t i = 0; i < A->cols(); ++i)
//      for (Eigen::SparseMatrix<double>::InnerIterator it(*A, i); it; ++it)
//        A_ub_->coeffRef(n_ub_ + it.row(), i) = it.value();

    b_ub_->conservativeResize(n_ub_ + n_add);
    b_ub_->tail(n_add) = *b;
  } else {
    A_ub_ = std::make_unique<Eigen::SparseMatrix<double>>(*A);
    b_ub_ = std::make_unique<Eigen::VectorXd>(*b);
  }

  ub_activity_.conservativeResize(n_ub_ + n_add);
  for (size_t i = 0; i < n_add; ++i)
    ub_activity_.coeffRef(n_ub_ + i) = false;
  ub_index_.conservativeResize(n_ub_ + n_add);

  n_ub_ += n_add;
}

void stuka::LP::LazyLinearProgram::removeConstr_ub(size_t index) {
  removeConstrs_ub(index, 1);
}

void stuka::LP::LazyLinearProgram::removeConstrs_ub(size_t index, size_t n_remove) {
  std::unique_ptr<Eigen::SparseMatrix<double>>
      A_ub = std::make_unique<Eigen::SparseMatrix<double>>(A_ub_->rows() - n_remove, n_dim_);
  A_ub->reserve(A_ub_->nonZeros());

  for (size_t i = 0; i < n_dim_; ++i) {
    A_ub->startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub_, i); it; ++it) {
      if (it.row() < index) A_ub->insertBack(it.row(), i) = it.value();
      else if (it.row() >= index + n_remove) A_ub->insertBack(it.row() - n_remove, i) = it.value();
    }
  }
  A_ub->finalize();

  std::unique_ptr<Eigen::VectorXd> b_ub = std::make_unique<Eigen::VectorXd>(b_ub_->size() - n_remove);
  b_ub->head(index) = b_ub_->head(index);
  b_ub->tail(b_ub->size() - index) = b_ub_->tail(b_ub_->size() - index - n_remove);

  b_ub_ = std::move(b_ub);

  size_t n_deactivated = 0;

  for (size_t i = 0; i < n_remove; ++i) {
    if (ub_activity_.coeff(index + i)) {
      prog_->removeConstr_ub(ub_index_.coeff(index + i)); // Remove constraint
      n_deactivated++;

      // Decrement indices
      for (size_t j = 0; j < ub_index_.size(); ++j)
        if (ub_activity_.coeff(j) && ub_index_.coeff(j) > ub_index_.coeff(index + i))
          ub_index_.coeffRef(j) -= 1;
    }
  }

  Eigen::Array<size_t, Eigen::Dynamic, 1> ub_index(n_ub_ - n_remove);
  ub_activity.head(index) = ub_activity_.head(index);
  ub_index.head(index) = ub_index_.head(index);
  ub_activity.tail(n_ub_ - index - n_remove) = ub_activity_.tail(n_ub_ - index - n_remove);
  ub_index.tail(ub_index_.size() - index - n_remove) = ub_index_.tail(ub_index_.size() - index - n_remove);
  ub_activity_ = ub_activity;
  ub_index_ = ub_index;

  n_active_ -= n_deactivated;
  n_ub_ -= n_remove;
}

void stuka::LP::LazyLinearProgram::addConstr_eq(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) {
  prog_->addConstr_eq(a, b);
}

void stuka::LP::LazyLinearProgram::addConstrs_eq(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                                                 const std::unique_ptr<Eigen::VectorXd> &b) {
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
  if (isActive(index)) return;

  std::unique_ptr<Eigen::VectorXd> a = std::make_unique<Eigen::VectorXd>(A_ub_->row(index));
  prog_->addConstr_ub(a, b_ub_->coeff(index));
  ub_activity_.coeffRef(index) = true;
  ub_index_.coeffRef(index) = n_active_;
  n_active_++;
}

void stuka::LP::LazyLinearProgram::setConstrsActive(size_t index, size_t n_active) {

  throw std::runtime_error("LP::LazyLinearProgram::setConstrsActive not implemented.");
  // TODO: Add check for whether constraints are already active.
  for (size_t i = 0; i < n_active; ++i)
    if (ub_activity_.coeff(index + i))
      throw std::runtime_error(
          "LP::LazyLinearProgram::setConstrsActive not implemented with existing active constraints.");

  std::unique_ptr<Eigen::SparseMatrix<double>>
      A = std::make_unique<Eigen::SparseMatrix<double>>(A_ub_->middleRows(index, n_active));
  std::unique_ptr<Eigen::VectorXd> b = std::make_unique<Eigen::VectorXd>(b_ub_->segment(index, n_active));
  prog_->addConstrs_ub(A, b);
  for (size_t i = 0; i < n_active; ++i) {
    ub_activity_.coeffRef(index + i) = true;
    ub_index_.coeffRef(index + i) = n_active_++;
  }
}

void stuka::LP::LazyLinearProgram::setConstrsActive(Eigen::Array<bool, Eigen::Dynamic, 1> activity) {
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

size_t stuka::LP::LazyLinearProgram::addViolatedConstraints(const Eigen::VectorXd &x) {
  if (!A_ub_ || x.size() == 0) return false;
  size_t added = 0;
  Eigen::Array<bool, Eigen::Dynamic, 1> violated = ((*A_ub_ * x - *b_ub_).array() > 0);
  for (size_t i = 0; i < violated.size(); ++i) {
    if (violated.coeff(i) && !ub_activity_.coeff(i)) {
      setConstrActive(i);
      added++;
  }
  }
  return added;
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
  Eigen::VectorXd dual = Eigen::VectorXd::Constant(n_ub_, 0);
  for (size_t i = 0; i < n_ub_; ++i)
    if (ub_activity_.coeff(i))
      dual.coeffRef(i) = active_dual.coeff(ub_index_.coeff(i));
  return dual;
}
