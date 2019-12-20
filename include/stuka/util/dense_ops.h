//
// Created by Avinash Madavan on 12/18/19.
//

#ifndef STUKA_UTIL_DENSE_OPS_H_
#define STUKA_UTIL_DENSE_OPS_H_

#include <Eigen/Core>

namespace stuka { namespace util { namespace DenseOps {

using ActiveSet = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
using MatRowPairD = std::pair<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
                              const Eigen::Matrix<bool, Eigen::Dynamic, 1>>;
using VecRowPairD = std::pair<const Eigen::VectorXd, const Eigen::Matrix<bool, Eigen::Dynamic, 1>>;

template<typename T, int Rows = Eigen::Dynamic, int Cols = Eigen::Dynamic>
Eigen::Matrix<T, Eigen::Dynamic, Cols> get_rows(const Eigen::Matrix<T, Rows, Cols> &mat,
                                                const Eigen::Matrix<bool, Rows, 1> &rows) {

  Eigen::Matrix<T, Eigen::Dynamic, Cols> ret(rows.count());

  size_t current_row = 0;
  for (size_t i = 0; i < mat.rows(); ++i)
    if (rows.coeff(i))
      ret.row(current_row++) = mat.row(i);

  return ret;
}

template<typename T, int Rows = Eigen::Dynamic, int Cols = Eigen::Dynamic>
Eigen::Matrix<T, Eigen::Dynamic, Cols> vstack(const std::vector<const Eigen::Matrix<T, Rows, Cols>> &mats) {

  size_t n_col = mats[0].cols(), n_row = 0;

  for (const Eigen::Matrix<T, Rows, Cols> &mat: mats)
    n_row += mat.rows();

  Eigen::Matrix<T, Eigen::Dynamic, Cols> ret(n_row, n_col);

  size_t row_offset = 0;
  for (const Eigen::Matrix<T, Rows, Cols> &mat: mats) {
    size_t rows = mat.rows();
    ret.rows(row_offset + rows) = mat;
    row_offset += rows;
  }

  return ret;
}

template<typename T, int Rows = Eigen::Dynamic, int Cols = Eigen::Dynamic>
Eigen::Matrix<T, Eigen::Dynamic, Cols> vstack_rows(const std::vector<std::pair<const Eigen::Matrix<T, Rows, Cols>,
                                                                               const Eigen::Matrix<bool, Rows, 1>>
> &mat_rows) {

  size_t n_mats = mat_rows.size(), n_col = mat_rows.begin()->first.cols(), n_row_out = 0;

  std::vector<Eigen::VectorXi> prev_row_counts(n_mats);
  size_t n = 0;

  for (const std::pair<const Eigen::Matrix<T, Rows, Cols>, const Eigen::Matrix<bool, Eigen::Dynamic, 1>>
        &mat_row:mat_rows) {
    const Eigen::Matrix<T, Rows, Cols> &mat = mat_row.first;
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
    n++;
  }

  Eigen::Matrix<T, Eigen::Dynamic, Cols> ret(n_row_out, n_col);

  for (size_t i = 0; i < n_col; ++i) {
    size_t row_offset = 0, n = 0;
    for (const std::pair<const Eigen::Matrix<T, Rows, Cols>, const Eigen::Matrix<bool, Eigen::Dynamic, 1>>
          &mat_row : mat_rows) {
      const Eigen::Matrix<T, Rows, Cols> &mat = mat_row.first;
      const Eigen::Matrix<bool, Eigen::Dynamic, 1> &rows = mat_row.second;

      for (size_t j = 0; j < mat.rows(); ++j)
        if (rows.coeff(j))
          ret.coeffRef(row_offset + prev_row_counts[n].coeff(j), i) = mat.coeff(j, i);

      n++;
      row_offset += rows.count();
    }
  }

  return ret;

}

}}}

#endif //STUKA_UTIL_DENSE_OPS_H_
