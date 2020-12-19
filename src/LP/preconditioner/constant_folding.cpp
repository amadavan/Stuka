//
// Created by Avinash Madavan on 2019-06-26.
//

#include <stuka/LP/preconditioner/constant_folding.h>

stuka::LP::ConstantFolding::ConstantFolding(const std::shared_ptr<stuka::LP::BaseLinearProgram> &next)
    : BasePreconditioner(next) {}

void stuka::LP::ConstantFolding::initialize(const stuka::LP::LinearProgram &prog) {
  n_dim_ = prog.c->size();
  n_ub_ = (prog.b_ub) ? prog.b_ub->size() : 0;
  n_eq_ = (prog.b_eq) ? prog.b_eq->size() : 0;

  A_ub_ = util::SparseOps::unique_copy(prog.A_ub);
  b_ub_ = util::DenseOps::unique_copy(prog.b_ub);
  A_eq_ = util::SparseOps::unique_copy(prog.A_eq);
  b_eq_ = util::DenseOps::unique_copy(prog.b_eq);

  is_constant_ = std::vector<bool>(n_dim_, false);
  constant_value_ = Eigen::VectorXd(n_dim_);
  if (prog.lb && prog.ub)
    for (size_t i = 0; i < n_dim_; ++i)
      if (prog.lb->coeff(i) == prog.ub->coeff(i)) {
        is_constant_[i] = true;
        constant_value_.coeffRef(i) = prog.lb->coeff(i);
      }

  n_constants_ = std::count(is_constant_.begin(), is_constant_.end(), true);

  LinearProgram lp_next;
  if (n_constants_ > 0) {
    size_t n_dim_next = n_dim_ - n_constants_;
    lp_next.c = std::make_unique<Eigen::VectorXd>(n_dim_next);
    if (prog.b_ub) {
      size_t n_ub = lp_next.b_ub->size();
      lp_next.A_ub = std::make_unique<Eigen::SparseMatrix<double>>(n_ub, n_dim_next);
      lp_next.A_ub->reserve(prog.A_ub->nonZeros());
      lp_next.b_ub = util::DenseOps::unique_copy(prog.b_ub);

      b_ub_shift_ = Eigen::VectorXd(prog.b_ub->size());
    }
    if (prog.b_eq) {
      size_t n_eq = lp_next.b_eq->size();
      lp_next.A_eq = std::make_unique<Eigen::SparseMatrix<double>>(n_eq, n_dim_next);
      lp_next.A_eq->reserve(prog.A_eq->nonZeros());
      lp_next.b_eq = util::DenseOps::unique_copy(prog.b_eq);

      b_eq_shift_ = Eigen::VectorXd(prog.b_eq->size());
    }
    lp_next.lb = std::make_unique<Eigen::VectorXd>(n_dim_next);
    lp_next.ub = std::make_unique<Eigen::VectorXd>(n_dim_next);

    size_t num_const = 0;
    for (size_t i = 0; i < n_dim_; ++i) {
      if (is_constant_[i]) {
        if (lp_next.b_ub) b_ub_shift_ -= lp_next.A_ub->col(i) * lp_next.lb->coeff(i);
        if (lp_next.b_eq) b_eq_shift_ -= lp_next.A_eq->col(i) * lp_next.lb->coeff(i);
        num_const++;
      } else {
        lp_next.c->coeffRef(i - num_const) = prog.c->coeff(i);
        if (prog.b_ub) {
          lp_next.A_ub->startVec(i - num_const);
          for (Eigen::SparseMatrix<double>::InnerIterator it_ub(*prog.A_ub, i); it_ub; ++it_ub)
            lp_next.A_ub->insertBack(it_ub.row(), i - num_const) = it_ub.value();
        }
        if (prog.b_eq) {
          lp_next.A_eq->startVec(i - num_const);
          for (Eigen::SparseMatrix<double>::InnerIterator it_eq(*prog.A_eq, i); it_eq; ++it_eq)
            lp_next.A_eq->insertBack(it_eq.row(), i - num_const) = it_eq.value();
        }
        lp_next.lb->coeffRef(i - num_const) = prog.lb->coeff(i);
        lp_next.ub->coeffRef(i - num_const) = prog.ub->coeff(i);
      }
    }

    if (lp_next.b_ub) {
      lp_next.A_ub->finalize();
      *lp_next.b_ub -= b_ub_shift_;
    }
    if (lp_next.b_eq) {
      lp_next.A_eq->finalize();
      *lp_next.b_eq -= b_eq_shift_;
    }
  } else {
    lp_next.c = util::DenseOps::unique_copy(prog.c);
    lp_next.A_ub = util::SparseOps::unique_copy(prog.A_ub);
    lp_next.b_ub = util::DenseOps::unique_copy(prog.b_ub);
    lp_next.A_eq = util::SparseOps::unique_copy(prog.A_eq);
    lp_next.b_eq = util::DenseOps::unique_copy(prog.b_eq);
    lp_next.lb = util::DenseOps::unique_copy(prog.lb);
    lp_next.ub = util::DenseOps::unique_copy(prog.ub);
  }

  next()->initialize(lp_next);
}

void stuka::LP::ConstantFolding::setObjective(const std::unique_ptr<Eigen::VectorXd> &c) {
  if (n_constants_ == 0) {
    next()->setObjective(c);
    return;
  }

  std::unique_ptr<Eigen::VectorXd> c_next = std::make_unique<Eigen::VectorXd>(n_dim_ - n_constants_);

  size_t num_const = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
    if (is_constant_[i]) num_const++;
    else c_next->coeffRef(i - num_const) = c->coeff(i);
  }

  next()->setObjective(c_next);
}

void stuka::LP::ConstantFolding::setRHS(const std::unique_ptr<Eigen::VectorXd> &b_ub,
                                        const std::unique_ptr<Eigen::VectorXd> &b_eq) {
  if (n_constants_ == 0) {
    next()->setRHS(b_ub, b_eq);
    return;
  }

  std::unique_ptr<Eigen::VectorXd> b_ub_ = std::make_unique<Eigen::VectorXd>(*b_ub);
  std::unique_ptr<Eigen::VectorXd> b_eq_ = std::make_unique<Eigen::VectorXd>(*b_eq);

  if (b_ub) *b_ub_ += b_ub_shift_;
  if (b_eq) *b_eq_ += b_eq_shift_;

  next()->setRHS(b_ub_, b_eq_);
}

void stuka::LP::ConstantFolding::setBounds(const std::unique_ptr<Eigen::VectorXd> &lb,
                                           const std::unique_ptr<Eigen::VectorXd> &ub) {
  throw std::runtime_error("LP::ConstantFolding::setBounds not implemented");
}

void stuka::LP::ConstantFolding::addVar(double c, const std::unique_ptr<Eigen::VectorXd> &a_ub,
                                        const std::unique_ptr<Eigen::VectorXd> &a_eq, double lb, double ub) {
  if (b_ub_) {
    A_ub_->conservativeResize(n_ub_, n_dim_ + 1);
    A_ub_->reserve(A_ub_->nonZeros() + n_dim_);
    for (size_t i = 0; i < n_ub_; ++i)
      if (abs(a_ub->coeff(i)) > 0) A_ub_->coeffRef(i, n_dim_) = a_ub->coeff(i);
  }

  if (b_eq_) {
    A_eq_->conservativeResize(n_eq_, n_dim_ + 1);
    A_eq_->reserve(A_eq_->nonZeros() + n_dim_);
    for (size_t i = 0; i < n_eq_; ++i)
      if (abs(a_eq->coeff(i)) > 0) A_eq_->coeffRef(i, n_dim_) = a_ub->coeff(i);
  }

  constant_value_.conservativeResize(n_dim_ + 1);
  constant_value_.coeffRef(n_dim_) = 0.;

  is_constant_.resize(n_dim_ + 1, false);

  if (lb == ub) {
    constant_value_.coeffRef(n_dim_) = lb;
    is_constant_[n_dim_] = true;
    n_constants_++;

    std::unique_ptr<Eigen::VectorXd> b_ub = std::make_unique<Eigen::VectorXd>(*b_ub_);
    std::unique_ptr<Eigen::VectorXd> b_eq = std::make_unique<Eigen::VectorXd>(*b_eq_);

    if (b_ub) {
      b_ub_shift_ -= *a_ub * lb;
      *b_ub += b_ub_shift_;
    }

    if (b_eq) {
      b_eq_shift_ -= *a_eq * lb;
      *b_eq += b_eq_shift_;
    }

    next()->setRHS(b_ub, b_eq);
  } else {
    next()->addVar(c, a_ub, a_eq, lb, ub);
  }

  n_dim_++;
}

void stuka::LP::ConstantFolding::addVars(const std::unique_ptr<Eigen::VectorXd> &c,
                                         const std::unique_ptr<Eigen::SparseMatrix<double>> &A_ub,
                                         const std::unique_ptr<Eigen::SparseMatrix<double>> &A_eq,
                                         const std::unique_ptr<Eigen::VectorXd> &lb,
                                         const std::unique_ptr<Eigen::VectorXd> &ub) {

  size_t n_add = (c) ? c->size() : (lb) ? lb->size() : (ub) ? ub->size() : (A_ub) ? A_ub->cols() : (A_eq) ? A_eq->cols()
                                                                                                          : 0;

  for (size_t i = 0; i < n_add; ++i) {
    std::unique_ptr<Eigen::VectorXd> a_ub = nullptr;
    std::unique_ptr<Eigen::VectorXd> a_eq = nullptr;

    if (b_ub_) {
      a_ub = std::make_unique<Eigen::VectorXd>(n_ub_);
      a_ub->setZero();
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub, i); it; ++it)
        a_ub->coeffRef(it.row()) = it.value();
    }

    if (b_eq_) {
      a_eq = std::make_unique<Eigen::VectorXd>(n_eq_);
      a_eq->setZero();
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_eq, i); it; ++it)
        a_ub->coeffRef(it.row()) = it.value();
    }

    addVar((c) ? c->coeff(i) : 0, a_ub, a_eq, (lb) ? lb->coeff(i) : -INF, (ub) ? ub->coeff(i) : INF);
  }
}

void stuka::LP::ConstantFolding::removeVar(size_t var) {
  if (b_ub_) {
    std::unique_ptr<Eigen::SparseMatrix<double>>
        A_ub = std::make_unique<Eigen::SparseMatrix<double>>(n_ub_, n_dim_ - 1);
    A_ub->reserve(A_ub_->nonZeros());
    for (size_t i = 0; i < n_dim_; ++i) {
      if (i == var) {
        if (is_constant_[i])
          b_ub_shift_ -= A_ub_->col(i) * constant_value_.coeff(i);
        continue;
      }
      size_t col = (i > var) ? i - 1 : i;
      A_ub->startVec(col);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub_, i); it; ++it)
        A_ub->insertBack(it.row(), col) = it.value();
      A_ub->finalize();
    }
    A_ub_ = std::move(A_ub);
  }

  if (b_eq_) {
    std::unique_ptr<Eigen::SparseMatrix<double>>
        A_eq = std::make_unique<Eigen::SparseMatrix<double>>(n_eq_, n_dim_ - 1);
    A_eq->reserve(A_eq_->nonZeros());
    for (size_t i = 0; i < n_dim_; ++i) {
      if (i == var) {
        if (is_constant_[i])
          b_eq_shift_ -= A_eq_->col(i) * constant_value_.coeff(i);
        continue;
      }
      size_t col = (i > var) ? i - 1 : i;
      A_eq->startVec(col);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_eq_, i); it; ++it)
        A_eq->insertBack(it.row(), col) = it.value();
      A_eq->finalize();
    }
    A_eq_ = std::move(A_eq);
  }

  if (is_constant_[var]) {
    std::unique_ptr<Eigen::VectorXd> b_ub = util::DenseOps::unique_copy(b_ub_);
    std::unique_ptr<Eigen::VectorXd> b_eq = util::DenseOps::unique_copy(b_eq_);

    if (b_ub_) *b_ub += b_ub_shift_;
    if (b_eq_) *b_eq += b_eq_shift_;

    next()->setRHS(b_ub, b_eq);
  } else {
    next()->removeVar(var);
  }

  is_constant_.erase(is_constant_.begin() + var);
  Eigen::VectorXd constant_value(n_dim_ - 1);
  constant_value.head(var) = constant_value_.head(var);
  constant_value.tail(n_dim_ - 1 - var) = constant_value_.tail(n_dim_ - var - 1);
  constant_value_ = constant_value;
  n_dim_--;
}

void stuka::LP::ConstantFolding::removeVars(size_t index, size_t n_remove) {
  for (size_t i = 0; i < n_remove; ++i)
    removeVar(index + i);
}

void stuka::LP::ConstantFolding::removeBackVars(size_t n_remove) {
  for (size_t i = 0; i < n_remove; ++i)
    removeVar(n_dim_ - 1 - n_remove + i);
}

void stuka::LP::ConstantFolding::addConstr_ub(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) {
  if (b_ub_) {
    A_ub_->conservativeResize(n_ub_ + 1, n_dim_);
    b_ub_->conservativeResize(n_ub_ + 1);
    b_ub_shift_.conservativeResize(n_ub_ + 1);
  } else {
    A_ub_ = std::make_unique<Eigen::SparseMatrix<double>>(1, n_dim_);
    b_ub_ = std::make_unique<Eigen::VectorXd>(1);
    b_ub_shift_ = Eigen::VectorXd(1);
  }

  for (size_t i = 0; i < n_dim_; ++i)
    if (abs(a->coeff(i)) > 1e-16)
      A_ub_->coeffRef(n_ub_, i) = a->coeff(i);
  b_ub_->coeffRef(n_ub_) = b;
  b_ub_shift_.coeffRef(n_ub_) = 0;

  if (n_constants_ == 0) {
    next()->addConstr_ub(a, b);
  } else {
    std::unique_ptr<Eigen::VectorXd> a_next = std::make_unique<Eigen::VectorXd>(n_dim_ - n_constants_);

    size_t num_const = 0;
    for (size_t i = 0; i < n_dim_; ++i) {
      if (is_constant_[i]) {
        b_ub_shift_.coeffRef(n_ub_) -= a->coeff(i) * constant_value_.coeff(i);
        num_const++;
      } else
        a_next->coeffRef(i - num_const) = A_ub_->coeff(n_ub_, i);
    }
    next()->addConstr_ub(a_next, b - b_ub_shift_.coeff(n_ub_));
  }

  n_ub_++;
}

void stuka::LP::ConstantFolding::addConstrs_ub(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                                               const std::unique_ptr<Eigen::VectorXd> &b) {
  size_t n_add = b->size();

  if (b_ub_) {
    A_ub_->conservativeResize(n_ub_ + n_add, n_dim_);
    b_ub_->conservativeResize(n_ub_ + n_add);
    b_ub_shift_.conservativeResize(n_ub_ + n_add);
  } else {
    A_ub_ = std::make_unique<Eigen::SparseMatrix<double>>(n_add, n_dim_);
    b_ub_ = std::make_unique<Eigen::VectorXd>(n_add);
    b_ub_shift_ = Eigen::VectorXd(n_add);
  }

  for (size_t i = 0; i < n_dim_; ++i)
    for (Eigen::SparseMatrix<double>::InnerIterator it_ub(*A, i); it_ub; ++it_ub)
      A_ub_->coeffRef(n_ub_ + it_ub.row(), i) = it_ub.value();
  b_ub_->tail(n_add) = *b;
  b_ub_shift_.tail(n_ub_).setZero();

  if (n_constants_ == 0) {
    next()->addConstrs_ub(A, b);
  } else {
    std::unique_ptr<Eigen::SparseMatrix<double>> A_next = std::make_unique<Eigen::SparseMatrix<double>>(n_add, n_dim_ -
        n_constants_);

    size_t num_const = 0;
    for (size_t i = 0; i < n_dim_; ++i) {
      if (is_constant_[i]) {
        b_ub_shift_.tail(n_add) -= A->col(i) * constant_value_.coeff(i);
        num_const++;
      } else {
        A_next->startVec(i - num_const);
        for (Eigen::SparseMatrix<double>::InnerIterator it_ub(*A, i); it_ub; ++it_ub)
          A_next->insertBack(it_ub.row(), i - num_const) = it_ub.value();
      }
    }
    A_next->finalize();
    std::unique_ptr<Eigen::VectorXd>
        b_next = std::make_unique<Eigen::VectorXd>(b_ub_->tail(n_add) - b_ub_shift_.tail(n_add));
    next()->addConstrs_ub(A_next, b_next);
  }

  n_ub_ += n_add;
}

void stuka::LP::ConstantFolding::removeConstr_ub(size_t index) {
  std::unique_ptr<Eigen::SparseMatrix<double>> A = std::make_unique<Eigen::SparseMatrix<double>>(n_ub_ - 1, n_dim_);
  std::unique_ptr<Eigen::VectorXd> b = std::make_unique<Eigen::VectorXd>(n_ub_ - 1);

  A->reserve(A_ub_->nonZeros());
  for (size_t i = 0; i < n_dim_; ++i) {
    A->startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it_ub(*A_ub_, i); it_ub; ++it_ub) {
      if (it_ub.row() == index) continue;
      A->insertBack(it_ub.row() - ((it_ub.row() > index) ? 1 : 0), i) = it_ub.value();
    }
  }
  A->finalize();

  b->head(index) = b_ub_->head(index);
  b->tail(n_dim_ - index - 1) = b_ub_->tail(n_dim_ - index - 1);

  A_ub_ = std::move(A);
  b_ub_ = std::move(b);

  b_ub_shift_.head(index) = b_ub_shift_.head(index);
  b_ub_shift_.segment(index, n_dim_ - index - 1) = b_ub_shift_.tail(n_dim_ - index - 1);
  b_ub_shift_.conservativeResize(n_ub_ - 1);

  next()->removeConstr_ub(index);

  n_ub_--;
}

void stuka::LP::ConstantFolding::removeConstrs_ub(size_t index, size_t n_remove) {
  std::unique_ptr<Eigen::SparseMatrix<double>> A = std::make_unique<Eigen::SparseMatrix<double>>(n_ub_ - n_remove,
                                                                                                 n_dim_);
  std::unique_ptr<Eigen::VectorXd> b = std::make_unique<Eigen::VectorXd>(n_ub_ - 1);

  A->reserve(A_ub_->nonZeros());
  for (size_t i = 0; i < n_dim_; ++i) {
    A->startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it_ub(*A_ub_, i); it_ub; ++it_ub) {
      if (it_ub.row() >= index && it_ub.row() < index + n_remove) continue;
      A->insertBack(it_ub.row() - ((it_ub.row() > index) ? n_remove : 0), i) = it_ub.value();
    }
  }
  A->finalize();

  b->head(index) = b_ub_->head(index);
  b->tail(n_dim_ - index - n_remove) = b_ub_->tail(n_dim_ - index - n_remove);

  A_ub_ = std::move(A);
  b_ub_ = std::move(b);

  b_ub_shift_.head(index) = b_ub_shift_.head(index);
  b_ub_shift_.segment(index + n_remove, n_dim_ - index - n_remove) = b_ub_shift_.tail(n_dim_ - index - n_remove);
  b_ub_shift_.conservativeResize(n_ub_ - n_remove);

  next()->removeConstrs_ub(index, n_remove);

  n_ub_ -= n_remove;
}

void stuka::LP::ConstantFolding::addConstr_eq(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) {
  if (b_eq_) {
    A_eq_->conservativeResize(n_eq_ + 1, n_dim_);
    b_eq_->conservativeResize(n_eq_ + 1);
    b_eq_shift_.conservativeResize(n_eq_ + 1);
  } else {
    A_eq_ = std::make_unique<Eigen::SparseMatrix<double>>(1, n_dim_);
    b_eq_ = std::make_unique<Eigen::VectorXd>(1);
    b_eq_shift_ = Eigen::VectorXd(1);
  }

  for (size_t i = 0; i < n_dim_; ++i)
    if (abs(a->coeff(i)) > 1e-16)
      A_eq_->coeffRef(n_eq_, i) = a->coeff(i);
  b_eq_->coeffRef(n_eq_) = b;
  b_eq_shift_.coeffRef(n_eq_) = 0;

  if (n_constants_ == 0) {
    next()->addConstr_eq(a, b);
  } else {
    std::unique_ptr<Eigen::VectorXd> a_next = std::make_unique<Eigen::VectorXd>(n_dim_ - n_constants_);

    size_t num_const = 0;
    for (size_t i = 0; i < n_dim_; ++i) {
      if (is_constant_[i]) {
        b_eq_shift_.coeffRef(n_eq_) -= a->coeff(i) * constant_value_.coeff(i);
        num_const++;
      } else
        a_next->coeffRef(i - num_const) = A_eq_->coeff(n_eq_, i);
    }
    next()->addConstr_eq(a_next, b - b_eq_shift_.coeff(n_eq_));
  }

  n_eq_++;
}

void stuka::LP::ConstantFolding::addConstrs_eq(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                                               const std::unique_ptr<Eigen::VectorXd> &b) {
  size_t n_add = b->size();

  if (b_eq_) {
    A_eq_->conservativeResize(n_eq_ + n_add, n_dim_);
    b_eq_->conservativeResize(n_eq_ + n_add);
    b_eq_shift_.conservativeResize(n_eq_ + n_add);
  } else {
    A_eq_ = std::make_unique<Eigen::SparseMatrix<double>>(n_add, n_dim_);
    b_eq_ = std::make_unique<Eigen::VectorXd>(n_add);
    b_eq_shift_ = Eigen::VectorXd(n_add);
  }

  for (size_t i = 0; i < n_dim_; ++i)
    for (Eigen::SparseMatrix<double>::InnerIterator it_eq(*A, i); it_eq; ++it_eq)
      A_eq_->coeffRef(n_eq_ + it_eq.row(), i) = it_eq.value();
  b_eq_->tail(n_add) = *b;
  b_eq_shift_.tail(n_eq_).setZero();

  if (n_constants_ == 0) {
    next()->addConstrs_eq(A, b);
  } else {
    std::unique_ptr<Eigen::SparseMatrix<double>> A_next = std::make_unique<Eigen::SparseMatrix<double>>(n_add, n_dim_ -
        n_constants_);

    size_t num_const = 0;
    for (size_t i = 0; i < n_dim_; ++i) {
      if (is_constant_[i]) {
        b_eq_shift_.tail(n_add) -= A->col(i) * constant_value_.coeff(i);
        num_const++;
      } else {
        A_next->startVec(i - num_const);
        for (Eigen::SparseMatrix<double>::InnerIterator it_eq(*A, i); it_eq; ++it_eq)
          A_next->insertBack(it_eq.row(), i - num_const) = it_eq.value();
      }
    }
    A_next->finalize();
    std::unique_ptr<Eigen::VectorXd>
        b_next = std::make_unique<Eigen::VectorXd>(b_eq_->tail(n_add) - b_eq_shift_.tail(n_add));
    next()->addConstrs_eq(A_next, b_next);
  }

  n_eq_ += n_add;
}

void stuka::LP::ConstantFolding::removeConstr_eq(size_t index) {
  std::unique_ptr<Eigen::SparseMatrix<double>> A = std::make_unique<Eigen::SparseMatrix<double>>(n_eq_ - 1, n_dim_);
  std::unique_ptr<Eigen::VectorXd> b = std::make_unique<Eigen::VectorXd>(n_eq_ - 1);

  A->reserve(A_eq_->nonZeros());
  for (size_t i = 0; i < n_dim_; ++i) {
    A->startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it_eq(*A_eq_, i); it_eq; ++it_eq) {
      if (it_eq.row() == index) continue;
      A->insertBack(it_eq.row() - ((it_eq.row() > index) ? 1 : 0), i) = it_eq.value();
    }
  }
  A->finalize();

  b->head(index) = b_eq_->head(index);
  b->tail(n_dim_ - index - 1) = b_eq_->tail(n_dim_ - index - 1);

  A_eq_ = std::move(A);
  b_eq_ = std::move(b);

  b_eq_shift_.head(index) = b_eq_shift_.head(index);
  b_eq_shift_.segment(index, n_dim_ - index - 1) = b_eq_shift_.tail(n_dim_ - index - 1);
  b_eq_shift_.conservativeResize(n_eq_ - 1);

  next()->removeConstr_eq(index);

  n_eq_--;
}

void stuka::LP::ConstantFolding::removeConstrs_eq(size_t index, size_t n_remove) {
  std::unique_ptr<Eigen::SparseMatrix<double>> A = std::make_unique<Eigen::SparseMatrix<double>>(n_eq_ - n_remove,
                                                                                                 n_dim_);
  std::unique_ptr<Eigen::VectorXd> b = std::make_unique<Eigen::VectorXd>(n_eq_ - 1);

  A->reserve(A_eq_->nonZeros());
  for (size_t i = 0; i < n_dim_; ++i) {
    A->startVec(i);
    for (Eigen::SparseMatrix<double>::InnerIterator it_eq(*A_eq_, i); it_eq; ++it_eq) {
      if (it_eq.row() >= index && it_eq.row() < index + n_remove) continue;
      A->insertBack(it_eq.row() - ((it_eq.row() > index) ? n_remove : 0), i) = it_eq.value();
    }
  }
  A->finalize();

  b->head(index) = b_eq_->head(index);
  b->tail(n_dim_ - index - n_remove) = b_eq_->tail(n_dim_ - index - n_remove);

  A_eq_ = std::move(A);
  b_eq_ = std::move(b);

  b_eq_shift_.head(index) = b_eq_shift_.head(index);
  b_eq_shift_.segment(index + n_remove, n_dim_ - index - n_remove) = b_eq_shift_.tail(n_dim_ - index - n_remove);
  b_eq_shift_.conservativeResize(n_eq_ - n_remove);

  next()->removeConstrs_eq(index, n_remove);

  n_eq_ -= n_remove;
}

Eigen::VectorXd stuka::LP::ConstantFolding::convertState(const Eigen::VectorXd &x) {
  Eigen::VectorXd x_(n_dim_ - n_constants_);

  size_t num_const = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
    if (is_constant_[i]) num_const++;
    else x_.coeffRef(i - num_const) = x.coeff(i);
  }

  return next()->convertState(x_);
}

Eigen::VectorXd stuka::LP::ConstantFolding::revertState(const Eigen::VectorXd &x) {
  Eigen::VectorXd xn = next()->revertState(x);
  Eigen::VectorXd x_(n_dim_);

  size_t num_const = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
    if (is_constant_[i]) {
      x_.coeffRef(i) = constant_value_.coeff(i);
      num_const++;
    } else {
      x_.coeffRef(i) = x.coeff(i - num_const);
    }
  }

  return x_;
}