//
// Created by Avinash Madavan on 2019-05-08.
//

#include <stuka/util/callback/save_hdf5.h>
#include <iostream>

stuka::util::callback::SaveHDF5::SaveHDF5(const std::string &filename, const size_t n_entries, const bool compress)
    : n_entries_(n_entries), compress_(compress) {
  file_ = H5::H5File(filename, H5F_ACC_TRUNC);
}

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

  H5::DSetCreatPropList plist_x, plist_ub, plist_eq;
  if (compress_) {
    unsigned szip_pixels_per_block = 16;
    hsize_t chunk_x[2] = {std::min((size_t) (125000 / x_dim_), n_entries_), x_dim_};
    plist_x.setChunk(2, chunk_x);
    plist_x.setSzip(H5_SZIP_NN_OPTION_MASK, szip_pixels_per_block);

    if (dual_ub_dim_ > 0) {
      hsize_t chunk_ub[2] = {std::min((size_t) (125000 / dual_ub_dim_), n_entries_), dual_ub_dim_};
      plist_ub.setChunk(2, chunk_ub);
      plist_ub.setSzip(H5_SZIP_NN_OPTION_MASK, szip_pixels_per_block);
    }

    if (dual_eq_dim_ > 0) {
      hsize_t chunk_eq[2] = {std::min((size_t) (125000 / dual_eq_dim_), n_entries_), dual_eq_dim_};
      plist_eq.setChunk(2, chunk_eq);
      plist_eq.setSzip(H5_SZIP_NN_OPTION_MASK, szip_pixels_per_block);
    }
  }

  ds_x_ = file_.createDataSet("x", H5::PredType::NATIVE_DOUBLE, x, plist_x);

  if (dual_ub_dim_ > 0)
    ds_dual_ub_ = file_.createDataSet("dual_ub", H5::PredType::NATIVE_DOUBLE, dual_ub, plist_ub);

  if (dual_eq_dim_ > 0)
    ds_dual_eq_ = file_.createDataSet("dual_eq", H5::PredType::NATIVE_DOUBLE, dual_eq, plist_eq);
}

void stuka::util::callback::SaveHDF5::callback(const stuka::OptimizeState state) {
  hsize_t offset[2] = {state.nit, 0};
  hsize_t count[2] = {1, x_dim_};
  hsize_t dim_mem[2] = {1, x_dim_};
  H5::DataSpace d_mem = H5::DataSpace(2, dim_mem);

  H5::DataSpace x = ds_x_.getSpace();
  x.selectHyperslab(H5S_SELECT_SET, count, offset, nullptr, nullptr);
  ds_x_.write(state.x.data(), H5::PredType::NATIVE_DOUBLE, d_mem, x);

  if (dual_ub_dim_ > 0) {
    count[1] = dual_ub_dim_;
    dim_mem[1] = dual_ub_dim_;
    d_mem = H5::DataSpace(2, dim_mem);

    H5::DataSpace dual_ub = ds_dual_ub_.getSpace();
    dual_ub.selectHyperslab(H5S_SELECT_SET, count, offset, nullptr, nullptr);
    ds_dual_ub_.write(state.dual_ub.data(), H5::PredType::NATIVE_DOUBLE, d_mem, dual_ub);
  }

  if (dual_eq_dim_ > 0) {
    count[1] = dual_eq_dim_;
    dim_mem[1] = dual_eq_dim_;
    d_mem = H5::DataSpace(2, dim_mem);

    H5::DataSpace dual_eq = ds_dual_eq_.getSpace();
    dual_eq.selectHyperslab(H5S_SELECT_SET, count, offset, nullptr, nullptr);
    ds_dual_eq_.write(state.dual_eq.data(), H5::PredType::NATIVE_DOUBLE, d_mem, dual_eq);
  }
}
