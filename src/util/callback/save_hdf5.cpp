//
// Created by Avinash Madavan on 2019-05-08.
//

#include <stuka/util/callback/save_hdf5.h>

stuka::util::callback::SaveHDF5::SaveHDF5(const std::string &filename, const size_t n_entries) :
    file_(H5::H5File(filename, H5F_ACC_TRUNC)), n_entries_(n_entries) {}

void stuka::util::callback::SaveHDF5::initialize(const stuka::OptimizeState state) {
  x_dim_ = (hsize_t) state.x.size();
  dual_ub_dim_ = (hsize_t) state.dual_ub.size();
  dual_eq_dim_ = (hsize_t) state.dual_eq.size();

  hsize_t x_dim[2] = {n_entries_, x_dim_};
  hsize_t dual_ub_dim[2] = {n_entries_, dual_ub_dim_};
  hsize_t dual_eq_dim[2] = {n_entries_, dual_eq_dim_};

  H5::DataSpace x = H5::DataSpace(2, x_dim);
  H5::DataSpace dual_ub = H5::DataSpace(2, dual_ub_dim);
  H5::DataSpace dual_eq = H5::DataSpace(2, dual_eq_dim);

  ds_x_ = file_.createDataSet("x", H5::PredType::NATIVE_DOUBLE, x);
  ds_dual_ub_ = file_.createDataSet("dual_ub", H5::PredType::NATIVE_DOUBLE, dual_ub);
  ds_dual_eq_ = file_.createDataSet("dual_eq", H5::PredType::NATIVE_DOUBLE, dual_eq);
}

void stuka::util::callback::SaveHDF5::callback(const stuka::OptimizeState state) {
  hsize_t offset[2] = {state.nit, 0};
  hsize_t count[2] = {1, x_dim_};
  hsize_t dim_mem[2] = {1, x_dim_};
  H5::DataSpace d_mem = H5::DataSpace(2, dim_mem);

  H5::DataSpace x = ds_x_.getSpace();
  x.selectHyperslab(H5S_SELECT_SET, count, offset, nullptr, nullptr);
  ds_x_.write(state.x.data(), H5::PredType::NATIVE_DOUBLE, d_mem, x);

  count[1] = dual_ub_dim_;
  dim_mem[1] = dual_ub_dim_;
  d_mem = H5::DataSpace(2, dim_mem);

  H5::DataSpace dual_ub = ds_x_.getSpace();
  dual_ub.selectHyperslab(H5S_SELECT_SET, count, offset, nullptr, nullptr);
  ds_dual_ub_.write(state.x.data(), H5::PredType::NATIVE_DOUBLE, d_mem, dual_ub);

  count[1] = dual_eq_dim_;
  dim_mem[1] = dual_eq_dim_;
  d_mem = H5::DataSpace(2, dim_mem);

  H5::DataSpace dual_eq = ds_x_.getSpace();
  dual_eq.selectHyperslab(H5S_SELECT_SET, count, offset, nullptr, nullptr);
  ds_dual_eq_.write(state.x.data(), H5::PredType::NATIVE_DOUBLE, d_mem, dual_eq);
}
