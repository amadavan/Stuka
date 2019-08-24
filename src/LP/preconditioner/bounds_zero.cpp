//
// Created by Avinash Madavan on 2019-07-06.
//

#include <stuka/LP/preconditioner/bounds_zero.h>
#include <iostream>

stuka::LP::BoundsZero::BoundsZero(const std::shared_ptr<stuka::LP::BaseLinearProgram> &next) : BasePreconditioner(
    next) {}


void stuka::LP::BoundsZero::initialize(const stuka::LP::LinearProgram &prog) {
  n_dim_ = prog.c->size();
  n_ub_ = (prog.b_ub) ? prog.b_ub->size() : 0;
  n_bounds_ = 0;
  n_split_ = 0;

  x_shift_ = Eigen::VectorXd(n_dim_);
  x_shift_.setZero();


  // Infer properties from bounds
  if (!prog.lb && !prog.ub) {
    n_split_ = n_dim_;
    x_type_ = std::vector<StateType>(n_dim_, StateType::UNBOUNDED);
  } else if (!prog.ub) {
    x_type_ = std::vector<StateType>(n_dim_, StateType::NORMAL);
    // No upper bound given (assumed infinity)
    for (size_t i = 0; i < n_dim_; ++i) {
      if (prog.lb->coeff(i) == -INF) {
        n_split_++;
        x_type_[i] = StateType::UNBOUNDED;
      } else {
        x_shift_.coeffRef(i) = prog.lb->coeff(i);
      }
    }
  } else if (!prog.lb) {
    x_type_ = std::vector<StateType>(n_dim_, StateType::NEGATIVE);
    // No lower bound given (assumed -infinity)
    for (size_t i = 0; i < n_dim_; ++i) {
      if (prog.ub->coeff(i) == INF) {
        n_split_++;
        x_type_[i] = StateType::UNBOUNDED;
      } else {
        x_shift_.coeffRef(i) = -prog.ub->coeff(i);
      }
    }
  } else {
    x_type_ = std::vector<StateType>(n_dim_, StateType::NORMAL);
    // Both upper and lower bounds are defined
    for (size_t i = 0; i < n_dim_; ++i) {
      if (prog.lb->coeff(i) == -INF && prog.ub->coeff(i) == INF) {
        // No bounds
        n_split_++;
        x_type_[i] = StateType::UNBOUNDED;
      } else if (prog.lb->coeff(i) != -INF && prog.ub->coeff(i) != INF) {
        // Compact
        n_bounds_++;
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
        throw std::runtime_error("LP::BoundsZero::initialize something is seriously wrong with the provided bounds.");
      }
    }
  }

  LinearProgram lp_next;

//  n_slack_ = n_ub_ + n_bounds_;
  size_t n_dim_tot = n_dim_ + n_split_;
  size_t n_eq = (prog.b_eq) ? prog.b_eq->size() : 0;

  lp_next.c = std::make_unique<Eigen::VectorXd>(n_dim_tot);
  lp_next.c->setZero();
  if (prog.b_eq) {
    lp_next.A_eq = std::make_unique<Eigen::SparseMatrix<double>>(n_eq, n_dim_tot);
    lp_next.A_eq->reserve(((prog.A_eq) ? prog.A_eq->nonZeros() : 0));
    lp_next.b_eq = std::make_unique<Eigen::VectorXd>(n_eq);
  }
  if (prog.b_ub || n_bounds_) {
    lp_next.A_ub = std::make_unique<Eigen::SparseMatrix<double>>(n_bounds_ + n_ub_, n_dim_tot);
    lp_next.A_ub->reserve(((prog.A_ub) ? prog.A_ub->nonZeros() : 0) + n_bounds_);
    lp_next.b_ub = std::make_unique<Eigen::VectorXd>(n_ub_ + n_bounds_);
  }
  lp_next.lb = std::make_unique<Eigen::VectorXd>(n_dim_tot);
  lp_next.lb->setZero();

  // Initialize constraint bound
  if (prog.b_eq) *lp_next.b_eq = *prog.b_eq - *prog.A_eq * x_shift_;
  if (prog.b_ub) lp_next.b_ub->tail(n_ub_) = *prog.b_ub - *prog.A_ub * x_shift_;

  size_t current_shift = 0, con_shift = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
    if (prog.b_eq) lp_next.A_eq->startVec(i + current_shift);
    if (prog.b_ub || n_bounds_) lp_next.A_ub->startVec(i + current_shift);
    switch (x_type_[i]) {
      case StateType::NORMAL:
        lp_next.c->coeffRef(i + current_shift) = prog.c->coeff(i);
        if (prog.b_eq)
          for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_eq, i); it; ++it)
            lp_next.A_eq->insertBack(it.row(), i + current_shift) = it.value();
        if (prog.b_ub || n_bounds_)
          for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_ub, i); it; ++it)
            lp_next.A_ub->insertBack(n_bounds_ + it.row(), i + current_shift) = it.value();
        break;
      case StateType::NEGATIVE:
        lp_next.c->coeffRef(i + current_shift) = -prog.c->coeff(i);
        if (prog.b_eq)
          for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_eq, i); it; ++it)
            lp_next.A_eq->insertBack(it.row(), i + current_shift) = -it.value();
        if (prog.b_ub || n_bounds_)
          for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_ub, i); it; ++it)
            lp_next.A_ub->insertBack(n_bounds_ + it.row(), i + current_shift) = -it.value();
        break;
      case StateType::COMPACT:
        lp_next.c->coeffRef(i + current_shift) = prog.c->coeff(i);

        // Add additional constraint from bound
        lp_next.A_ub->insertBack(con_shift, i + current_shift) = 1.;
        lp_next.b_ub->coeffRef(con_shift) = prog.ub->coeff(i);
        con_shift++;

        // Nominal constraints
        if (prog.b_eq)
          for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_eq, i); it; ++it)
            lp_next.A_eq->insertBack(it.row(), i + current_shift) = it.value();
        if (prog.b_ub || n_bounds_)
          for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_ub, i); it; ++it)
            lp_next.A_ub->insertBack(n_bounds_ + it.row(), i + current_shift) = it.value();
        break;
      case StateType::UNBOUNDED:
        lp_next.c->coeffRef(i + current_shift) = prog.c->coeff(i);
        if (prog.b_eq)
          for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_eq, i); it; ++it)
            lp_next.A_eq->insertBack(it.row(), i + current_shift) = it.value();
        if (prog.b_ub || n_bounds_)
          for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_ub, i); it; ++it)
            lp_next.A_ub->insertBack(n_bounds_ + it.row(), i + current_shift) = it.value();

        current_shift++;

        if (prog.b_eq) lp_next.A_eq->startVec(i + current_shift);
        if (prog.b_ub || n_bounds_) lp_next.A_ub->startVec(i + current_shift);

        lp_next.c->coeffRef(i + current_shift) = -prog.c->coeff(i);
        if (prog.b_eq)
          for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_eq, i); it; ++it)
            lp_next.A_eq->insertBack(it.row(), i + current_shift) = -it.value();
        if (prog.b_ub || n_bounds_)
          for (Eigen::SparseMatrix<double>::InnerIterator it(*prog.A_ub, i); it; ++it)
            lp_next.A_ub->insertBack(n_bounds_ + it.row(), i + current_shift) = -it.value();
        break;
    }
  }

  lp_next.A_eq->finalize();
  lp_next.A_ub->finalize();

  next()->initialize(lp_next);
}

void stuka::LP::BoundsZero::setObjective(const std::shared_ptr<Eigen::VectorXd> &c) {

}

void stuka::LP::BoundsZero::setRHS(const std::shared_ptr<Eigen::VectorXd> &b_ub,
                                   const std::shared_ptr<Eigen::VectorXd> &b_eq) {

}

void stuka::LP::BoundsZero::setBounds(const std::shared_ptr<Eigen::VectorXd> &lb,
                                      const std::shared_ptr<Eigen::VectorXd> &ub) {

}

void
stuka::LP::BoundsZero::addVar(double c, std::shared_ptr<Eigen::VectorXd> a_ub, std::shared_ptr<Eigen::VectorXd> a_eq,
                              double lb, double ub) {

}

void
stuka::LP::BoundsZero::addVars(std::shared_ptr<Eigen::VectorXd> c, std::shared_ptr<Eigen::SparseMatrix<double>> A_ub,
                               std::shared_ptr<Eigen::SparseMatrix<double>> A_eq, std::shared_ptr<Eigen::VectorXd> lb,
                               std::shared_ptr<Eigen::VectorXd> ub) {

}

void stuka::LP::BoundsZero::removeVar(size_t var) {

}

void stuka::LP::BoundsZero::removeVars(size_t index, size_t n_remove) {

}

void stuka::LP::BoundsZero::removeBackVars(size_t n_remove) {

}

void stuka::LP::BoundsZero::addConstr_ub(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) {

}

void stuka::LP::BoundsZero::addConstrs_ub(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                                          const std::shared_ptr<Eigen::VectorXd> &b) {

}

void stuka::LP::BoundsZero::removeConstr_ub(size_t index) {

}

void stuka::LP::BoundsZero::removeConstrs_ub(size_t index, size_t n_remove) {

}

void stuka::LP::BoundsZero::addConstr_eq(const std::shared_ptr<Eigen::VectorXd> &a, const double &b) {

}

void stuka::LP::BoundsZero::addConstrs_eq(const std::shared_ptr<Eigen::SparseMatrix<double>> &A,
                                          const std::shared_ptr<Eigen::VectorXd> &b) {

}

void stuka::LP::BoundsZero::removeConstr_eq(size_t index) {

}

void stuka::LP::BoundsZero::removeConstrs_eq(size_t index, size_t n_remove) {

}

Eigen::VectorXd stuka::LP::BoundsZero::convertState(const Eigen::VectorXd &x) {
  Eigen::VectorXd x_slack(n_dim_ + n_split_);
  x_slack.setZero();

  size_t current_shift = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
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

    }
  }

  return next()->convertState(x_slack);
}

Eigen::VectorXd stuka::LP::BoundsZero::revertState(const Eigen::VectorXd &x) {
  Eigen::VectorXd xn = next()->revertState(x);
  Eigen::VectorXd x_ = x_shift_;

  size_t current_shift = 0;
  for (size_t i = 0; i < n_dim_; ++i) {
    switch (x_type_[i]) {
      case StateType::COMPACT:
      case StateType::NORMAL:
        x_.coeffRef(i) += xn.coeff(i + current_shift);
        break;
      case StateType::NEGATIVE :
        x_.coeffRef(i) -= xn.coeff(i + current_shift);
        break;
      case StateType::UNBOUNDED :
        x_.coeffRef(i) += xn.coeff(i + current_shift) - xn.coeff(i + current_shift + 1);
        current_shift += 1;
        break;
    }
  }

  return x_;
}