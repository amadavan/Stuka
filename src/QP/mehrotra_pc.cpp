#include <stuka/QP/mehrotra_pc.h>

// TODO: implement quadratic costs

stuka::QP::MehrotraPC::MehrotraPC(const LP::LinearProgram &lp, const Options &opts) : BaseLPSolver(lp, opts), prog_(lp), opts_(opts) {
  // Set number of variables in original problem
  n_base_ = lp.c->size();
  n_compact_ = 0;
  n_free_ = 0;

  m_eq_ = (lp.b_eq != nullptr) ? lp.b_eq->size() : 0;
  m_ub_ = (lp.b_ub != nullptr) ? lp.b_ub->size() : 0;

  // Evaluate variable types
  if (lp.lb != nullptr && lp.ub != nullptr) {
    for (int i = 0; i < n_base_; ++i) {
      if (lp.lb->coeff(i) != -INF && lp.ub->coeff(i) != INF) {
        variable_types_[i] = COMPACT;
        n_compact_++;
      }
      else if (lp.lb->coeff(i) != -INF)
        variable_types_[i] = (lp.lb->coeff(i) == 0) ? STANDARD : LOWER;
      else if (lp.ub->coeff(i) != INF)
        variable_types_[i] = UPPER;
      else {
        variable_types_[i] = FREE;
        n_free_++;
      }
    }
  } else if (lp.lb != nullptr) {
    for (int i = 0; i < n_base_; ++i)
      variable_types_[i] = (lp.lb->coeff(i) == 0) ? STANDARD : LOWER;
  } else if (lp.ub != nullptr) {
    for (int i = 0; i < n_base_; ++i)
      variable_types_[i] = UPPER;
  } else { // lp.lb == lp.ub == nullptr
    for (int i = 0; i < n_base_; ++i)
      variable_types_[i] = FREE;
    n_free_ = n_base_;
  }

  // Initialize matrices and vectors
  n_ = n_base_ + n_free_ + m_ub_;
  A_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(m_eq_ + m_ub_, n_);
  b_ = Eigen::VectorXd(m_eq_ + m_ub_);
  c_ = Eigen::VectorXd(n_);
  u_ = Eigen::VectorXd(n_compact_);
  bound_ = Eigen::VectorXd(n_base_);

  // Populate matrices and vectors
  int nnz = ((m_eq_ > 0) ? lp.A_eq->nonZeros() : 0) + ((m_ub_ > 0) ? lp.A_ub->nonZeros() : 0);
  double density = nnz / (double) (n_base_ * (m_eq_ + m_ub_));
  // A_.reserve(nnz + m_ub_ + (int) ceil(density * n_free_ * (m_eq_ + m_ub_)));
  if (m_eq_ > 0) b_.head(m_eq_) = *lp.b_eq;
  if (m_ub_ > 0) b_.tail(m_ub_) = *lp.b_ub;

  std::map<int, int> index_offset; // Map original index to free-var-adjusted

  int ind_free = 0;     ///< Number of free variables already added
  int ind_compact = 0;  ///< Number of compact variables already added
  for (auto it = variable_types_.cbegin(); it != variable_types_.cend(); ++it) {
    index_offset[it->first] = it->first + ind_free;
    if (it->second == STANDARD) {
      c_.coeffRef(index_offset[it->first]) = lp.c->coeff(it->first);
    } else if (it->second == LOWER) {
      c_.coeffRef(index_offset[it->first]) = lp.c->coeff(it->first);
      b_.head(m_eq_) -= lp.A_eq->col(it->first) * lp.lb->coeff(it->first);
      b_.tail(m_ub_) -= lp.A_ub->col(it->first) * lp.lb->coeff(it->first);
      bound_.coeffRef(it->first) = lp.lb->coeff(it->first);
    } else if (it->second == UPPER) {
      c_.coeffRef(index_offset[it->first]) = -lp.c->coeff(it->first);
      b_.head(m_eq_) -= lp.A_eq->col(it->first) * lp.ub->coeff(it->first);
      b_.tail(m_ub_) -= lp.A_ub->col(it->first) * lp.ub->coeff(it->first);
      bound_.coeffRef(it->first) = lp.ub->coeff(it->first);
    } else if (it->second == COMPACT) {
      c_.coeffRef(index_offset[it->first]) = lp.c->coeff(it->first);
      b_.head(m_eq_) -= lp.A_eq->col(it->first) * lp.lb->coeff(it->first);
      b_.tail(m_ub_) -= lp.A_ub->col(it->first) * lp.lb->coeff(it->first);
      bound_.coeffRef(it->first) = lp.lb->coeff(it->first);
      u_.coeffRef(ind_compact) = lp.ub->coeff(it->first) - lp.lb->coeff(it->first);
      z_index_[ind_compact] = index_offset[it->first];
      ind_compact++;
    } else if (it->second == FREE) {
      c_.coeffRef(index_offset[it->first]) = lp.c->coeff(it->first);
      c_.coeffRef(index_offset[it->first] + 1) = -lp.c->coeff(it->first);
      ind_free++;
    } else { // This should never be reached
    }
  }
  
  // Allocate A using triplet insertion
  std::vector<Eigen::Triplet<double>> A_triplets;
  A_triplets.reserve(nnz + m_ub_ + (int) ceil(density * n_free_ * (m_eq_ + m_ub_)));
  for (int i = 0; i < n_base_; ++i) {
    if (m_eq_ > 0) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(*lp.A_eq, i); it; ++it) {
        A_triplets.push_back({it.row(), index_offset[it.col()], it.value()});
        if (variable_types_[i] == FREE)
          A_triplets.push_back({it.row(), index_offset[it.col()] + 1, -it.value()});
      }
    }
    if (m_ub_ > 0) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(*lp.A_ub, i); it; ++it) {
        A_triplets.push_back({m_eq_ + it.row(), index_offset[it.col()], it.value()});
        if (variable_types_[i] == FREE)
          A_triplets.push_back({m_eq_ + it.row(), index_offset[it.col()] + 1, -it.value()});
      }
    }
  }

  for (int i = 0; i < m_ub_; ++i)
    A_triplets.push_back({m_eq_ + i, n_base_ + n_free_ + i, 1});

  A_.setFromTriplets(A_triplets.begin(), A_triplets.end());
  A_.makeCompressed();

  c_norm_ = c_.norm() + 1;
  b_norm_ = b_.norm() + 1;
  u_norm_ = u_.norm() + 1;

  computeInitialPoint();
}

void stuka::QP::MehrotraPC::iterate() {

  mu_ = x_.dot(s_) / n_;
  if (n_compact_ > 0) mu_ += z_.dot(w_) / n_compact_;

  rhs_.Rc = c_ - A_.transpose() * y_ - s_;
  rhs_.Rb = b_ - A_ * x_;
  rhs_.Ru = u_ - z_;
  rhs_.Rxs = -x_.cwiseProduct(s_);
  rhs_.Rzw = -z_.cwiseProduct(w_);
  
  // Adjustments for compact variables
  for (auto it = z_index_.cbegin(); it != z_index_.cend(); ++it) {
    rhs_.Rc.coeffRef(it->second) += w_.coeff(it->first);
    rhs_.Ru.coeffRef(it->first) -= x_.coeff(it->second);
  }

  // Compute the affine scaling step
  Step step_aff = solveKKT(rhs_, true);
  StepSize stepsize_aff = computeStepSize(step_aff, 1.);

  // Set scaling parameters
  double mu_aff = (x_ + stepsize_aff.primal * step_aff.dx).dot(s_ + stepsize_aff.dual * step_aff.ds) / n_;
  if (n_compact_ > 0) mu_aff += (z_ + stepsize_aff.primal * step_aff.dz).dot(w_ + stepsize_aff.dual * step_aff.dw) / n_compact_;
  double sigma = pow(mu_aff / mu_, 3);

  // Adjust the rhs according to the affine step
  rhs_.Rxs -= step_aff.dx.cwiseProduct(step_aff.ds);
  rhs_.Rxs.array() += sigma * mu_;

  rhs_.Rzw -= step_aff.dz.cwiseProduct(step_aff.dw);
  rhs_.Rzw.array() += sigma * mu_;

  // Compute the center-corrector step
  Step step_cc = solveKKT(rhs_, false);
  StepSize stepsize_cc = computeStepSize(step_cc, 0.99995);

  // Update the variables
  x_ += stepsize_cc.primal * step_cc.dx;
  z_ += stepsize_cc.primal * step_cc.dz;
  y_ += stepsize_cc.dual * step_cc.dy;
  s_ += stepsize_cc.dual * step_cc.ds;
  w_ += stepsize_cc.dual * step_cc.dw;
}

bool stuka::QP::MehrotraPC::terminate() {
  error_ = rhs_.Rc.norm() / c_norm_ + rhs_.Rb.norm() / b_norm_ + x_.dot(s_) / n_;
  if (n_compact_ > 0) error_ += rhs_.Ru.norm() / u_norm_ + z_.dot(w_) / n_compact_;
  return error_ < opts_.tol;
}

const stuka::OptimizeState stuka::QP::MehrotraPC::getState() {
  OptimizeState state;
  state.fun = c_.dot(x_);
  state.x = Eigen::VectorXd(n_base_);
  state.dual_eq = y_.head(m_eq_);
  state.dual_ub = y_.tail(m_ub_);
  state.dual_x_lb = Eigen::VectorXd::Zero(n_base_);
  state.dual_x_ub = Eigen::VectorXd::Zero(n_base_);
  state.status;

  // Adjust state to original format
  size_t n_free = 0;
  size_t n_compact = 0;
  for (auto it = variable_types_.cbegin(); it != variable_types_.cend(); ++it) {
    if (it->second == STANDARD) {
      state.x.coeffRef(it->first) = x_.coeff(it->first + n_free);
      state.dual_x_lb.coeffRef(it->first) = s_.coeff(it->first + n_free);
    } else if (it->second == LOWER) {
      state.x.coeffRef(it->first) = x_.coeff(it->first + n_free) + bound_.coeff(it->first);
      state.dual_x_lb.coeffRef(it->first) = s_.coeff(it->first + n_free);
    } else if (it->second == UPPER) {
      state.x.coeffRef(it->first) = bound_.coeff(it->first) - x_.coeff(it->first + n_free);
      state.dual_x_ub.coeffRef(it->first) = s_.coeff(it->first + n_free);
    } else if (it->second == FREE) {
      state.x.coeffRef(it->first) = x_.coeff(it->first + n_free) - x_.coeff(it->first + n_free + 1);
      n_free++;
    } else if (it->second == COMPACT) {
      state.x.coeffRef(it->first) = x_.coeff(it->first + n_free);
      state.dual_x_lb.coeffRef(it->first) = s_.coeff(it->first + n_free);
      state.dual_x_ub.coeffRef(it->first) = w_.coeff(n_compact);
      n_compact++;
    } else {
      // TODO: should be unreachable
    }
  }

  return state;
}

void stuka::QP::MehrotraPC::computeInitialPoint() {
  // Compute factorization
  Eigen::SparseMatrix<double, Eigen::RowMajor> AA_t = A_ * A_.transpose();
  ADA_t.analyzePattern(AA_t);
  ADA_t.factorize(AA_t);

  // Select x as min norm(x) s.t. Ax = b
  x_ = A_.transpose() * ADA_t.solve(b_);
  
  // Select z as z = u - x
  z_ = u_;
  for (auto it = z_index_.cbegin(); it != z_index_.cend(); ++it)
    z_.coeffRef(it->first) -= x_.coeff(it->second);

  // Select y, s, w as min norm(s = w) s.t. A.T y + s - w = c
  y_ = ADA_t.solve((Eigen::SparseMatrix<double, Eigen::ColMajor>) A_) * c_;
  s_ = c_ - A_.transpose() * y_;
  w_ = Eigen::VectorXd::Zero(n_compact_);
  for (auto it = z_index_.cbegin(); it != z_index_.cend(); ++it) {
    if (s_.coeff(it->second) != 0) {
      w_.coeffRef(it->first) =  s_.coeff(it->second);
      s_.coeffRef(it->second) = 0;
    }
  }

  double dx = std::max(-1.5 * x_.minCoeff(), 0.);
  double dz = (n_compact_ > 0) ? std::max(-1.5 * z_.minCoeff(), 0.) : 0;
  double ds = std::max(-1.5 * s_.minCoeff(), 0.);
  double dw = (n_compact_ > 0) ? std::max(-1.5 * w_.minCoeff(), 0.) : 0;

  double pdct = 0.5 * (x_.array() + dx).matrix().dot((s_.array() + ds).matrix());
  double dx_c = dx + pdct / (s_.sum() + n_ * ds);
  double ds_c = ds + pdct / (x_.sum() + n_ * dx);

  double pdct_c = 0.5 * (z_.array() + dz).matrix().dot((w_.array() + dw).matrix());
  double dz_c = dz + pdct_c / (w_.sum() + n_ * dw);
  double dw_c = dw + pdct_c / (z_.sum() + n_ * dz);

  // Adjust initial point
  x_.array() += dx_c;
  z_.array() += dz_c;
  s_.array() += ds_c;
  w_.array() += dw_c;
}

stuka::QP::MehrotraPC::Step stuka::QP::MehrotraPC::solveKKT(stuka::QP::MehrotraPC::RHS rhs, bool resolve) {
  // Compute theta^-1 = [X^-1 S + Z^-1 W]^-1
  Eigen::VectorXd theta = s_.cwiseQuotient(x_);

  Eigen::VectorXd z_inv_w = w_.cwiseQuotient(z_);
  for (auto it = z_index_.cbegin(); it != z_index_.cend(); ++it)
    theta.coeffRef(it->second) += z_inv_w.coeff(it->first);

  Eigen::VectorXd theta_inv = theta.cwiseInverse();

  // Solve the linear system
  if (resolve) ADA_t.factorize(A_ * theta_inv.asDiagonal() * A_.transpose());

  // Compute the step
  Step step = {};

  Eigen::VectorXd common_terms = rhs.Rc - rhs.Rxs.cwiseQuotient(x_);
  if (n_compact_ > 0) common_terms += (rhs.Rzw - w_.cwiseProduct(rhs.Ru)).cwiseQuotient(z_);
  step.dy = ADA_t.solve(rhs.Rb + A_ * theta_inv.asDiagonal() * common_terms);
  step.dx = theta_inv.asDiagonal() * (A_.transpose() * step.dy - common_terms);
  step.ds = (rhs.Rxs - s_.cwiseProduct(step.dx)).cwiseQuotient(x_);
  if (n_compact_ > 0) {
    step.dz = rhs.Ru;
    for (auto it = z_index_.cbegin(); it != z_index_.cend(); ++it)
      step.dz.coeffRef(it->first) -= step.dx.coeff(it->second);
    step.dw = (rhs.Rzw - w_.cwiseProduct(step.dz)).cwiseQuotient(z_);
  }

  return step;
}

stuka::QP::MehrotraPC::StepSize stuka::QP::MehrotraPC::computeStepSize(stuka::QP::MehrotraPC::Step step, double scale) {
  StepSize step_size = {1., 1.};

  for (int i = 0; i < n_; ++i) {
    if (step.dx.coeff(i) < 0) {
      double step_size_i = scale * -x_.coeff(i)/step.dx.coeff(i);
      step_size.primal = std::min(step_size_i, step_size.primal);
    }
    if (step.ds.coeff(i) < 0) {
      double step_size_i = scale * -s_.coeff(i)/step.ds.coeff(i);
      step_size.dual = std::min(step_size_i, step_size.dual);
    }
  }

  for (int i = 0; i < n_compact_; ++i) {
    if (step.dz.coeff(i) < 0) {
      double step_size_i = scale * -z_.coeff(i)/step.dz.coeff(i);
      step_size.primal = std::min(step_size_i, step_size.primal);
    }
    if (step.dw.coeff(i) < 0) {
      double step_size_i = scale * -w_.coeff(i)/step.dw.coeff(i);
      step_size.dual = std::min(step_size_i, step_size.dual);
    }
  }

  return step_size;
}

stuka::LP::BaseLinearProgram &stuka::QP::MehrotraPC::getLP() {
  throw std::runtime_error("QP::MehrotraPC currently does not support modification of LP");
}