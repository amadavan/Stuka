//
// Created by Avinash Madavan on 2019-03-02.
//

#include <stuka/LP/slack_lp.h>

#include <iostream>

void stuka::LP::SlackLinearProgram::initialize(const stuka::LP::LinearProgram &prog) {
  long n_dim = prog.c->size();
  long n_con_ub = (prog.b_ub) ? prog.b_ub->size() : 0;
  long n_con_eq = (prog.b_eq) ? prog.b_eq->size() : 0;

  x_shift_ = Eigen::VectorXd(n_dim);
  x_shift_.setZero();

  n_split_ = 0;
  n_const_ = 0;
  n_add_ = 0;

  // Infer properties from bounds
  if (!prog.lb && !prog.ub) {
    n_split_ = n_dim;
    x_type_ = std::vector<StateType>(n_dim, StateType::UNBOUNDED);
  } else if (!prog.ub) {
    x_type_ = std::vector<StateType>(n_dim, StateType::NORMAL);
    // No upper bound given (assumed infinity)
    for (size_t i = 0; i < n_dim; ++i) {
      if (prog.lb->coeff(i) == -INF) {
        n_split_++;
        x_type_[i] = StateType::UNBOUNDED;
      } else {
        x_shift_.coeffRef(i) = prog.lb->coeff(i);
      }
    }
  } else if (!prog.lb) {
    x_type_ = std::vector<StateType>(n_dim, StateType::NEGATIVE);
    // No lower bound given (assumed -infinity)
    for (size_t i = 0; i < n_dim; ++i) {
      if (prog.ub->coeff(i) == INF) {
        n_split_++;
        x_type_[i] = StateType::UNBOUNDED;
      } else {
        x_shift_.coeffRef(i) = -prog.ub->coeff(i);
      }
    }
  } else {
    x_type_ = std::vector<StateType>(n_dim, StateType::NORMAL);
    // Both upper and lower bounds are defined
    for (size_t i = 0; i < n_dim; ++i) {
      if (prog.lb->coeff(i) == -INF && prog.ub->coeff(i) == INF) {
        // No bounds
        n_split_++;
        x_type_[i] = StateType::UNBOUNDED;
      } else if (prog.lb->coeff(i) == prog.ub->coeff(i)) {
        // Constant
        n_const_++;
        x_type_[i] = StateType::CONSTANT;
      } else if (prog.lb->coeff(i) != -INF && prog.ub->coeff(i) != INF) {
        // Compact
        n_add_++;
        x_type_[i] = StateType::COMPACT;
        x_shift_.coeffRef(i) = prog.lb->coeff(i);
      } else if (prog.lb->coeff(i) != -INF) {
        // No upper bound
        x_shift_.coeffRef(i) = prog.lb->coeff(i);
      } else if (prog.ub->coeff(i) != INF) {
        // No lower bound
        x_type_[i] = StateType::NEGATIVE;
        x_shift_.coeffRef(i) = -prog.ub->coeff(i);
      } else {
        // This code should never be reached (the above cases should be exhaustive)
        throw std::runtime_error("LP::SlackLinearProgram something is seriously wrong with the provided bounds.");
      }
    }
  }

  // Construct problem
  n_slack_ = n_con_ub + n_add_;
  n_dim_ = n_dim + n_split_ - n_const_ + n_slack_;
  n_con_ = n_con_eq + n_slack_;

  c_ = std::make_unique<Eigen::VectorXd>(n_dim_);
  c_->setZero();
  A_ = std::make_unique<Eigen::SparseMatrix<double>>(n_con_, n_dim_);
  b_ = std::make_unique<Eigen::VectorXd>(n_con_);

  // Initialize constraint bound
  b_->segment(n_add_, n_con_eq) = *prog.b_eq;
  b_->tail(n_con_ub) = *prog.b_ub;

  // Compute the effect of the shift in x for constraint bound
  b_->segment(n_add_, n_con_eq) -= *prog.A_eq * x_shift_;
  b_->tail(n_con_ub) -= *prog.A_ub * x_shift_;

  // Construct A and c
  A_->reserve(prog.A_ub->nonZeros() + prog.A_eq->nonZeros() + n_add_);
  size_t current_shift = 0, con_shift = 0;

  for (size_t i = 0; i < n_dim; ++i) {
    if (x_type_[i] != StateType::CONSTANT) A_->startVec(i + current_shift);
    switch (x_type_[i]) {
      case StateType::NORMAL:
        c_->coeffRef(i + current_shift) = prog.c->coeff(i);
        for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_eq, i); it; ++it)
          A_->insertBack(n_add_ + it.row(), i + current_shift) = it.value();
        for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_ub, i); it; ++it)
          A_->insertBack(n_add_ + n_con_eq + it.row(), i + current_shift) = it.value();
        break;
      case StateType::NEGATIVE:
        c_->coeffRef(i + current_shift) = -prog.c->coeff(i);
        for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_eq, i); it; ++it)
          A_->insertBack(n_add_ + it.row(), i + current_shift) = -it.value();
        for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_ub, i); it; ++it)
          A_->insertBack(n_add_ + n_con_eq + it.row(), i + current_shift) = -it.value();
        break;
      case StateType::COMPACT:
        c_->coeffRef(i + current_shift) = prog.c->coeff(i);

        // Add additional constraint from bound
        A_->insertBack(con_shift, i + current_shift) = 1.;
        b_->coeffRef(con_shift) = prog.ub->coeff(i);
        con_shift++;

        // Nominal constraints
        for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_eq, i); it; ++it)
          A_->insertBack(n_add_ + it.row(), i + current_shift) = it.value();
        for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_ub, i); it; ++it)
          A_->insertBack(n_add_ + n_con_eq + it.row(), i + current_shift) = it.value();
        break;
      case StateType::UNBOUNDED:
        c_->coeffRef(i + current_shift) = prog.c->coeff(i);
        for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_eq, i); it; ++it)
          A_->insertBack(n_add_ + it.row(), i + current_shift) = it.value();
        for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_ub, i); it; ++it)
          A_->insertBack(n_add_ + n_con_eq + it.row(), i + current_shift) = it.value();

        current_shift++;

        A_->startVec(i + current_shift);
        c_->coeffRef(i + current_shift) = -prog.c->coeff(i);
        for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_eq, i); it; ++it)
          A_->insertBack(n_add_ + it.row(), i + current_shift) = -it.value();
        for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_ub, i); it; ++it)
          A_->insertBack(n_add_ + n_con_eq + it.row(), i + current_shift) = -it.value();
        break;
      case StateType::CONSTANT:
        current_shift -= 1;
        break;
    }
  }

  // Add slack variables to A_ for additional variables from bounds
  for (size_t i = 0; i < n_add_; ++i) {
    A_->startVec(n_dim + n_split_ - n_const_ + i);
    A_->insertBack(i, n_dim + n_split_ - n_const_ + i) = 1.;
  }

  // Add slack variables to A_ for constraints
  for (size_t i = 0; i < n_con_ub; ++i) {
    A_->startVec(n_dim + n_split_ - n_const_ + n_add_ + i);
    A_->insertBack(n_add_ + n_con_eq + i, n_dim + n_split_ - n_const_ + n_add_ + i) = 1.;
  }

  A_->finalize();
}

void stuka::LP::SlackLinearProgram::setObjective(const std::shared_ptr<Eigen::VectorXd> &c) {
  throw std::runtime_error("LP::SlackLinearProgram::setObjective not implemented");

}

void stuka::LP::SlackLinearProgram::setRHS(const std::shared_ptr<Eigen::VectorXd> &b_ub,
                                           const std::shared_ptr<Eigen::VectorXd> &b_eq) {
  throw std::runtime_error("LP::SlackLinearProgram::setRHS not implemented");
}

void stuka::LP::SlackLinearProgram::setBounds(const std::shared_ptr<Eigen::VectorXd> &lb,
                                              const std::shared_ptr<Eigen::VectorXd> &ub) {
  throw std::runtime_error("LP::SlackLinearProgram::setBounds not implemented");
}

void stuka::LP::SlackLinearProgram::addVar(double c, std::shared_ptr<Eigen::VectorXd> a_ub,
                                           std::shared_ptr<Eigen::VectorXd> a_eq, double lb, double ub) {
  throw std::runtime_error("LP::SlackLinearProgram::addVar not implemented");
}

void stuka::LP::SlackLinearProgram::addVars(std::shared_ptr<Eigen::VectorXd> c,
                                            std::shared_ptr<Eigen::SparseMatrix<double>> A_ub,
                                            std::shared_ptr<Eigen::SparseMatrix<double>> A_eq,
                                            std::shared_ptr<Eigen::VectorXd> lb, std::shared_ptr<Eigen::VectorXd> ub) {
  throw std::runtime_error("LP::SlackLinearProgram::addVars not implemented");
}

void stuka::LP::SlackLinearProgram::removeVar(size_t var) {
  throw std::runtime_error("LP::SlackLinearProgram::removeVar not implemented");
}

void stuka::LP::SlackLinearProgram::removeVars(size_t index, size_t n_remove) {
  throw std::runtime_error("LP::SlackLinearProgram::removeVars not implemented");
}

void stuka::LP::SlackLinearProgram::removeBackVars(size_t n_remove) {
  throw std::runtime_error("LP::SlackLinearProgram::removeBackVars not implemented");
}

void stuka::LP::SlackLinearProgram::addConstr_ub(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) {
  throw std::runtime_error("LP::SlackLinearProgram::addConstr_ub not implemented");
}

void stuka::LP::SlackLinearProgram::addConstrs_ub(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                                                  const std::shared_ptr<Eigen::VectorXd> &b) {
  throw std::runtime_error("LP::SlackLinearProgram::addConstrs_ub not implemented");
}

void stuka::LP::SlackLinearProgram::removeConstr_ub(size_t index) {
  throw std::runtime_error("LP::SlackLinearProgram::removeConstr_ub not implemented");
}

void stuka::LP::SlackLinearProgram::removeConstrs_ub(size_t index, size_t n_remove) {
  throw std::runtime_error("LP::SlackLinearProgram::removeConstrs_ub not implemented");
}

void stuka::LP::SlackLinearProgram::addConstr_eq(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) {
  throw std::runtime_error("LP::SlackLinearProgram::addConstr_eq not implemented");
}

void stuka::LP::SlackLinearProgram::addConstrs_eq(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                                                  const std::shared_ptr<Eigen::VectorXd> &b) {
  throw std::runtime_error("LP::SlackLinearProgram::addConstrs_eq not implemented");
}

void stuka::LP::SlackLinearProgram::removeConstr_eq(size_t index) {
  throw std::runtime_error("LP::SlackLinearProgram::removeConstr_eq not implemented");
}

void stuka::LP::SlackLinearProgram::removeConstrs_eq(size_t index, size_t n_remove) {
  throw std::runtime_error("LP::SlackLinearProgram::removeConstrs_eq not implemented");
}

Eigen::VectorXd stuka::LP::SlackLinearProgram::convertState(const Eigen::VectorXd &x) {
  Eigen::VectorXd x_slack(n_dim_);
  x_slack.setZero();

  size_t current_shift = 0, con_shift = 0;
  for (size_t i = 0; i < n_dim_ - n_split_ + n_const_ - n_slack_; ++i) {
    switch (x_type_[i]) {
      case StateType::COMPACT:
      case StateType::NORMAL:
        x_slack.coeffRef(i + current_shift) = x.coeff(i) + x_shift_.coeff(i);
        break;
      case StateType::NEGATIVE:
        x_slack.coeffRef(i + current_shift) = -x.coeff(i) + x_shift_.coeff(i);
        break;
      case StateType::UNBOUNDED:
        if (x.coeff(i) > 0)
          x_slack.coeffRef(i + current_shift) = x.coeff(i) + x_shift_.coeff(i);
        else
          x_slack.coeffRef(i + current_shift + 1) = -x.coeff(i) + x_shift_.coeff(i);
        current_shift++;
        break;
      case StateType::CONSTANT:
        current_shift--;
        break;
    }
  }

  return x_slack;
}

Eigen::VectorXd stuka::LP::SlackLinearProgram::revertState(const Eigen::VectorXd &x_slack) {
  Eigen::VectorXd x(x_shift_);

  size_t current_shift = 0;
  for (size_t i = 0; i < n_dim_ - n_split_ + n_const_ - n_slack_; ++i) {
    switch (x_type_[i]) {
      case StateType::COMPACT:
      case StateType::NORMAL:
        x.coeffRef(i) += x_slack.coeff(i + current_shift);
        break;
      case StateType::NEGATIVE :
        x.coeffRef(i) -= x_slack.coeff(i + current_shift);
        break;
      case StateType::CONSTANT :
        current_shift -= 1;
        break;
      case StateType::UNBOUNDED :
        x.coeffRef(i) += x_slack.coeff(i + current_shift) - x_slack[i + current_shift + 1];
        current_shift += 1;
        break;
    }
  }

  return x;
}
