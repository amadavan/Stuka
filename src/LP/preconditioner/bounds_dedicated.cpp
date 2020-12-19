//
// Created by Avinash Madavan on 2019-07-03.
//

#include <stuka/LP/preconditioner/bounds_dedicated.h>

stuka::LP::BoundsDedicated::BoundsDedicated(const std::shared_ptr<stuka::LP::BaseLinearProgram> &next)
    : BasePreconditioner(next) {}

void stuka::LP::BoundsDedicated::initialize(const stuka::LP::LinearProgram &prog) {
  n_dim_ = prog.c->size();
  n_ub_ = (prog.b_ub) ? prog.b_ub->size() : 0;
  n_bounds_ = 0;

  std::vector<stateType> x_type;

  // Infer properties from bounds
  if (!prog.lb && !prog.ub) {
    x_type = std::vector<stateType>(n_dim_, stateType::UNBOUNDED);
  } else if (!prog.ub) {
    x_type = std::vector<stateType>(n_dim_, stateType::NORMAL);
    // No upper bound given (assumed infinity)
    for (size_t i = 0; i < n_dim_; ++i) {
      if (prog.lb->coeff(i) == -INF)
        x_type[i] = stateType::UNBOUNDED;
      else
        n_bounds_++;
    }
  } else if (!prog.lb) {
    x_type = std::vector<stateType>(n_dim_, stateType::NEGATIVE);
    // No lower bound given (assumed -infinity)
    for (size_t i = 0; i < n_dim_; ++i) {
      if (prog.ub->coeff(i) == INF)
        x_type[i] = stateType::UNBOUNDED;
      else
        n_bounds_++;
    }
  } else {
    x_type = std::vector<stateType>(n_dim_, stateType::NORMAL);
    // Both upper and lower bounds are defined
    for (size_t i = 0; i < n_dim_; ++i) {
      if (prog.lb->coeff(i) == -INF && prog.ub->coeff(i) == INF) {
        // No bounds
        x_type[i] = stateType::UNBOUNDED;
      } else if (prog.lb->coeff(i) != -INF && prog.ub->coeff(i) != INF) {
        // Compact
        x_type[i] = stateType::COMPACT;
        n_bounds_ += 2;
      } else if (prog.lb->coeff(i) != -INF) {
        // No upper bound
        n_bounds_ += 1;
      } else if (prog.ub->coeff(i) != INF) {
        // No lower bound
        x_type[i] = stateType::NEGATIVE;
        n_bounds_ += 1;
      } else {
        // This code should never be reached (the above cases should be exhaustive)
        throw std::runtime_error("LP::SlackLinearProgram something is seriously wrong with the provided bounds.");
      }
    }
  }

  A_ub_ = util::SparseOps::unique_copy(prog.A_ub);
  b_ub_ = util::DenseOps::unique_copy(prog.b_ub);

  A_bounds_ = std::make_unique<Eigen::SparseMatrix<double>>(n_bounds_, n_dim_);
  b_bounds_ = std::make_unique<Eigen::VectorXd>(n_bounds_);

  size_t n_bound = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
    A_bounds_->startVec(i);
    switch (x_type[i]) {
      case stateType::UNBOUNDED:break;
      case stateType::NORMAL:A_bounds_->insertBack(n_bound, i) = -1.;
        b_bounds_->coeffRef(n_bound) = -prog.lb->coeffRef(i);
        n_bound++;
        break;
      case stateType::NEGATIVE:A_bounds_->insertBack(n_bound, i) = 1.;
        b_bounds_->coeffRef(n_bound) = prog.ub->coeffRef(i);
        n_bound++;
        break;
      case stateType::COMPACT:A_bounds_->insertBack(n_bound, i) = -1.;
        b_bounds_->coeffRef(n_bound) = -prog.lb->coeffRef(i);
        n_bound++;
        A_bounds_->insertBack(n_bound, i) = 1.;
        b_bounds_->coeffRef(n_bound) = prog.ub->coeffRef(i);
        n_bound++;
        break;
    }
  }

  LinearProgram lp_next;
  lp_next.c = util::DenseOps::unique_copy(prog.c);
  if (n_bounds_ > 0) {
    lp_next.A_ub = std::make_unique<Eigen::SparseMatrix<double>>(n_bounds_ + n_ub_, n_dim_);
    lp_next.b_ub = std::make_unique<Eigen::VectorXd>(n_bounds_ + n_ub_);
    lp_next.b_ub->head(n_bounds_) = *b_bounds_;
    lp_next.b_ub->tail(n_ub_) = *b_ub_;

    for (size_t i = 0; i < n_dim_; ++i) {
      lp_next.A_ub->startVec(i);
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_bounds_, i); it; ++it)
        lp_next.A_ub->insertBack(it.row(), i) = it.value();
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub_, i); it; ++it)
        lp_next.A_ub->insertBack(n_bounds_ + it.row(), i) = it.value();
    }
    lp_next.A_ub->finalize();

    lp_next.A_eq = util::SparseOps::unique_copy(prog.A_eq);
    lp_next.b_eq = util::DenseOps::unique_copy(prog.b_eq);
  } else {
    lp_next.c = util::DenseOps::unique_copy(prog.c);
    lp_next.A_ub = util::SparseOps::unique_copy(prog.A_ub);
    lp_next.b_ub = util::DenseOps::unique_copy(prog.b_ub);
    lp_next.A_eq = util::SparseOps::unique_copy(prog.A_eq);
    lp_next.b_eq = util::DenseOps::unique_copy(prog.b_eq);
  }

  next()->initialize(lp_next);
}

void stuka::LP::BoundsDedicated::setObjective(const std::unique_ptr<Eigen::VectorXd> &c) {
  next()->setObjective(c);
}

void stuka::LP::BoundsDedicated::setRHS(const std::unique_ptr<Eigen::VectorXd> &b_ub,
                                        const std::unique_ptr<Eigen::VectorXd> &b_eq) {
  b_ub_ = util::DenseOps::unique_copy(b_ub);

  std::unique_ptr<Eigen::VectorXd> b_ub_next = std::make_unique<Eigen::VectorXd>(n_bounds_ + n_ub_ + 1);
  b_ub_next->head(n_bounds_);
  b_ub_next->tail(n_ub_) = *b_ub;

  next()->setRHS(b_ub_next, b_eq);

  n_ub_++;
}

void stuka::LP::BoundsDedicated::setBounds(const std::unique_ptr<Eigen::VectorXd> &lb,
                                           const std::unique_ptr<Eigen::VectorXd> &ub) {
  throw std::runtime_error("LP::BoundsDedicated::setBounds not implemented");
}

void stuka::LP::BoundsDedicated::addVar(double c, const std::unique_ptr<Eigen::VectorXd> &a_ub,
                                        const std::unique_ptr<Eigen::VectorXd> &a_eq, double lb, double ub) {
  if (lb == INF && ub == INF) {
    next()->addVar(c, a_ub, a_eq, lb, ub);
  } else {
    next()->addVar(c, a_ub, a_eq, lb, ub);

  }

  n_dim_++;

  // TODO: add next
}

// TODO: addVars needs to add constraints for bounds.
void stuka::LP::BoundsDedicated::addVars(const std::unique_ptr<Eigen::VectorXd> &c,
                                         const std::unique_ptr<Eigen::SparseMatrix<double>> &A_ub,
                                         const std::unique_ptr<Eigen::SparseMatrix<double>> &A_eq,
                                         const std::unique_ptr<Eigen::VectorXd> &lb,
                                         const std::unique_ptr<Eigen::VectorXd> &ub) {

  size_t
      n_add = (c) ? c->size() : (lb) ? lb->size() : (ub) ? ub->size() : (A_ub) ? A_ub->cols() : (A_eq) ? A_eq->cols()
                                                                                                       : 0;

  if (b_ub_) {
    A_ub_->conservativeResize(n_ub_, n_dim_ + n_add);
    for (size_t i = 0; i < n_add; ++i)
      for (Eigen::SparseMatrix<double>::InnerIterator it(*A_ub, i); it; ++it)
        A_ub_->coeffRef(it.row(), n_dim_ + i) = it.value();
  }

  n_dim_ += n_add;
}

void stuka::LP::BoundsDedicated::removeVar(size_t var) {

}

void stuka::LP::BoundsDedicated::removeVars(size_t index, size_t n_remove) {

}

void stuka::LP::BoundsDedicated::removeBackVars(size_t n_remove) {

}

void stuka::LP::BoundsDedicated::addConstr_ub(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) {

}

void stuka::LP::BoundsDedicated::addConstrs_ub(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                                               const std::unique_ptr<Eigen::VectorXd> &b) {

}

void stuka::LP::BoundsDedicated::removeConstr_ub(size_t index) {

}

void stuka::LP::BoundsDedicated::removeConstrs_ub(size_t index, size_t n_remove) {

}

void stuka::LP::BoundsDedicated::addConstr_eq(const std::unique_ptr<Eigen::VectorXd> &a, const double &b) {
  next()->addConstr_eq(a, b);
}

void stuka::LP::BoundsDedicated::addConstrs_eq(const std::unique_ptr<Eigen::SparseMatrix<double>> &A,
                                               const std::unique_ptr<Eigen::VectorXd> &b) {
  next()->addConstrs_eq(A, b);
}

void stuka::LP::BoundsDedicated::removeConstr_eq(size_t index) {
  next()->removeConstr_eq(index);
}

void stuka::LP::BoundsDedicated::removeConstrs_eq(size_t index, size_t n_remove) {
  next()->removeConstrs_eq(index, n_remove);
}

Eigen::VectorXd stuka::LP::BoundsDedicated::convertState(const Eigen::VectorXd &x) {
  return next()->convertState(x);
}

Eigen::VectorXd stuka::LP::BoundsDedicated::revertState(const Eigen::VectorXd &x) {
  return next()->revertState(x);
}
