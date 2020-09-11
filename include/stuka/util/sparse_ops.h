//
// Created by Avinash Madavan on 12/14/19.
//

#ifndef STUKA_UTIL_SPARSE_OPS_H_
#define STUKA_UTIL_SPARSE_OPS_H_

#include <Eigen/SparseCore>

namespace stuka { namespace util { namespace SparseOps {

using MatRowPairD = std::pair<const Eigen::SparseMatrix<double>, const Eigen::Matrix<bool, Eigen::Dynamic, 1>>;

template<typename T>
Eigen::SparseMatrix<T> get_rows(const Eigen::SparseMatrix<T> mat, const Eigen::Matrix<bool, Eigen::Dynamic, 1> rows) {
  size_t n_row_in = mat.rows(), n_col = mat.cols(), n_row_out = rows.count();

  Eigen::VectorXi prev_row_count(n_row_in);

  int count_row = 0;
  for (size_t i = 0; i < n_row_in; ++i) {
    prev_row_count.coeffRef(i) = count_row;
    if (rows.coeff(i)) ++count_row;
  }

  Eigen::SparseMatrix<T> ret(n_row_out, n_col);
  ret.reserve(mat.nonZeros());

  for (size_t i = 0; i < n_col; ++i) {
    ret.startVec(i);
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, i); it; ++it)
      if (rows.coeff(it.row()))
        ret.insertBack(prev_row_count.coeff(it.row()), i) = it.value();
  }

  ret.finalize();

  return ret;
}

template<typename T>
Eigen::SparseMatrix<T> vstack(std::initializer_list<const Eigen::SparseMatrix<T>> mats) {

  // Iterate once to get system properties
  size_t n_row = 0, n_col = mats.begin()->cols(), nnz = 0;
  for (const Eigen::SparseMatrix<T> &mat : mats) {
    n_row += mat.rows();
    nnz += mat.nonZeros();
  }

  Eigen::SparseMatrix<T> ret(n_row, n_col);
  ret.reserve(nnz);

  // Iterate to generate the full matrix
  for (size_t i = 0; i < n_col; ++i) {
    ret.startVec(i);
    size_t row_offset = 0;
    for (const Eigen::SparseMatrix<T> &mat : mats) {
      for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, i); it; ++it)
        ret.insertBack(row_offset + it.row(), i) = it.value();
      row_offset += mat.rows();
    }
  }

  ret.finalize();

  return ret;
}

template<typename T>
Eigen::SparseMatrix<T> vstack_rows(const std::vector<
std::pair<const Eigen::SparseMatrix<T>, const Eigen::Matrix<bool, Eigen::Dynamic, 1>>> &mat_rows) {

  size_t n_mats = mat_rows.size(), n_col = mat_rows.begin()->first.cols(), n_row_out = 0, nnz = 0;

  std::vector<Eigen::VectorXi> prev_row_counts(n_mats);
  size_t n = 0;

  for (const std::pair<const Eigen::SparseMatrix<T>, const Eigen::Matrix<bool, Eigen::Dynamic, 1>>
        &mat_row : mat_rows) {
    const Eigen::SparseMatrix<T> &mat = mat_row.first;
    const Eigen::Matrix<bool, Eigen::Dynamic, 1> &rows = mat_row.second;

    size_t n_row_in = mat.rows();

    Eigen::VectorXi prev_row_count(n_row_in);

    int count_row = 0;
    for (size_t i = 0; i < n_row_in; ++i) {
      prev_row_count.coeffRef(i) = count_row;
      if (rows.coeff(i)) ++count_row;
    }

    prev_row_counts[n] = prev_row_count;

    n_row_out += rows.count();
    nnz += mat.nonZeros();
    n++;
  }

  Eigen::SparseMatrix<T> ret(n_row_out, n_col);
  ret.reserve(nnz);

  for (size_t i = 0; i < n_col; ++i) {
    ret.startVec(i);

    size_t row_offset = 0, n = 0;
    for (const std::pair<const Eigen::SparseMatrix<T>, const Eigen::Matrix<bool, Eigen::Dynamic, 1>>
          &mat_row : mat_rows) {
      const Eigen::SparseMatrix<T> &mat = mat_row.first;
      const Eigen::Matrix<bool, Eigen::Dynamic, 1> &rows = mat_row.second;

      for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, i); it; ++it)
        if (rows.coeff(it.row()))
          ret.insertBack(row_offset + prev_row_counts[n].coeff(it.row()), i) = it.value();

      n++;
      row_offset += rows.count();
    }
  }

  ret.finalize();

  return ret;
}

}}}

#endif //STUKA_UTIL_SPARSE_OPS_H_
