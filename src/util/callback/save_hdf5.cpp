//
// Created by Avinash Madavan on 2019-05-08.
//

#include <stuka/util/callback/save_hdf5.h>
#include <iostream>

stuka::util::callback::SaveHDF5::SaveHDF5(const std::string &filename, const size_t n_iter, const bool compress)
    : compress_(compress) {
  file_ = H5Fcreate(filename.data(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hid_t fprops = H5Fget_access_plist(file_);
  int mdc;
  size_t ccelems;
  size_t ccnbytes;
  double w0;
  H5Pget_cache(fprops, &mdc, &ccelems, &ccnbytes, &w0);
  H5Pset_cache(fprops, 0, 512, HDF5_CHUNK_CACHE * sizeof(double) * 3, w0);
  H5Pclose(fprops);

  n_p_ = ((size_t) (n_iter / 10000) > 0) ? (n_iter / 10000) : 1;
  n_entries_ = ((size_t) (n_iter / 10000) > 0) ? 10001 : n_iter + 1;
  index_ = 0;
}

stuka::util::callback::SaveHDF5::~SaveHDF5() {

  H5Dclose(ds_nit_);
  H5Dclose(ds_x_);
  if (dual_ub_dim_ > 0) H5Dclose(ds_dual_ub_);
  if (dual_eq_dim_ > 0) H5Dclose(ds_dual_eq_);
  H5Fclose(file_);

}

void stuka::util::callback::SaveHDF5::initialize(const stuka::OptimizeState state) {
  x_dim_ = (size_t) state.x.size();
  dual_ub_dim_ = (size_t) state.dual_ub.size();
  dual_eq_dim_ = (size_t) state.dual_eq.size();

  hsize_t nit_dim[1] = {n_entries_};
  hsize_t x_dim[2] = {n_entries_, x_dim_};
  hsize_t dual_ub_dim[2] = {n_entries_, dual_ub_dim_};
  hsize_t dual_eq_dim[2] = {n_entries_, dual_eq_dim_};

  hid_t nit = H5Screate_simple(1, nit_dim, nullptr);
  hid_t x = H5Screate_simple(2, x_dim, nullptr);
  hid_t dual_ub = H5Screate_simple(2, dual_ub_dim, nullptr);
  hid_t dual_eq = H5Screate_simple(2, dual_eq_dim, nullptr);

  hid_t plist_x = H5Pcreate(H5P_DATASET_CREATE);
  hid_t plist_ub = H5Pcreate(H5P_DATASET_CREATE);
  hid_t plist_eq = H5Pcreate(H5P_DATASET_CREATE);

  if (compress_) {
    size_t szip_pixels_per_block = (n_entries_ < 16) ? 2 : 16;

    hsize_t chunk_x[2] = {std::min((size_t) (HDF5_CHUNK_CACHE / x_dim_), n_entries_), x_dim_};
    H5Pset_chunk(plist_x, 2, chunk_x);
    H5Pset_szip(plist_x, H5_SZIP_NN_OPTION_MASK, szip_pixels_per_block);

    if (dual_ub_dim_ > 0) {
      hsize_t chunk_ub[2] = {std::min((size_t) (HDF5_CHUNK_CACHE / dual_ub_dim_), n_entries_), dual_ub_dim_};
      H5Pset_chunk(plist_ub, 2, chunk_ub);
      H5Pset_szip(plist_ub, H5_SZIP_NN_OPTION_MASK, szip_pixels_per_block);
    }

    if (dual_eq_dim_ > 0) {
      hsize_t chunk_eq[2] = {std::min((size_t) (HDF5_CHUNK_CACHE / dual_eq_dim_), n_entries_), dual_eq_dim_};
      H5Pset_chunk(plist_eq, 2, chunk_eq);
      H5Pset_szip(plist_eq, H5_SZIP_NN_OPTION_MASK, szip_pixels_per_block);
    }
  }

  ds_nit_ = H5Dcreate(file_, "nit", H5T_NATIVE_ULONG, nit, H5P_LINK_CREATE_DEFAULT, H5P_DATASET_CREATE_DEFAULT,
                      H5P_DATASET_ACCESS_DEFAULT);
  ds_x_ = H5Dcreate(file_, "x", H5T_NATIVE_DOUBLE, x, H5P_LINK_CREATE_DEFAULT, plist_x, H5P_DATASET_ACCESS_DEFAULT);
  if (dual_ub_dim_ > 0)
    ds_dual_ub_ = H5Dcreate(file_, "dual_ub", H5T_NATIVE_DOUBLE, dual_ub, H5P_LINK_CREATE_DEFAULT, plist_ub,
                            H5P_DATASET_ACCESS_DEFAULT);
  if (dual_eq_dim_ > 0)
    ds_dual_eq_ = H5Dcreate(file_, "dual_eq", H5T_NATIVE_DOUBLE, dual_eq, H5P_LINK_CREATE_DEFAULT, plist_eq,
                            H5P_DATASET_ACCESS_DEFAULT);

  H5Pclose(plist_x);
  H5Pclose(plist_ub);
  H5Pclose(plist_eq);

  H5Sclose(nit);
  H5Sclose(x);
  H5Sclose(dual_ub);
  H5Sclose(dual_eq);
}

void stuka::util::callback::SaveHDF5::callback(const stuka::OptimizeState state) {
  if (state.nit % n_p_ == 0)
    save(state);
}

void stuka::util::callback::SaveHDF5::finish(const stuka::OptimizeState state) {
  if (n_p_ > 1)
    save(state);
}

void stuka::util::callback::SaveHDF5::save(const stuka::OptimizeState state) {

  hsize_t offset_nit[1] = {index_};
  hsize_t count_nit[1] = {1};
  hsize_t dim_mem_nit[1] = {1};

  hid_t d_mem = H5Screate_simple(1, dim_mem_nit, nullptr);
  hid_t nit = H5Dget_space(ds_nit_);
  H5Sselect_hyperslab(nit, H5S_SELECT_SET, offset_nit, nullptr, count_nit, nullptr);
  H5Dwrite(ds_nit_, H5T_NATIVE_ULONG, d_mem, nit, H5P_DATASET_XFER_DEFAULT, &state.nit);

  H5Sclose(d_mem);
  H5Sclose(nit);

  hsize_t offset[2] = {index_, 0};
  hsize_t count[2] = {1, x_dim_};
  hsize_t dim_mem[2] = {1, x_dim_};

  d_mem = H5Screate_simple(2, dim_mem, nullptr);
  hid_t x = H5Dget_space(ds_x_);
  H5Sselect_hyperslab(x, H5S_SELECT_SET, offset, nullptr, count, nullptr);
  H5Dwrite(ds_x_, H5T_NATIVE_DOUBLE, d_mem, x, H5P_DATASET_XFER_DEFAULT, state.x.data());

  H5Sclose(d_mem);
  H5Sclose(x);

  if (dual_ub_dim_ > 0) {
    count[1] = dual_ub_dim_;
    dim_mem[1] = dual_ub_dim_;

    d_mem = H5Screate_simple(2, dim_mem, nullptr);
    hid_t dual_ub = H5Dget_space(ds_dual_ub_);
    H5Sselect_hyperslab(dual_ub, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    H5Dwrite(ds_dual_ub_, H5T_NATIVE_DOUBLE, d_mem, dual_ub, H5P_DATASET_XFER_DEFAULT, state.dual_ub.data());

    H5Sclose(d_mem);
    H5Sclose(dual_ub);
  }

  if (dual_eq_dim_ > 0) {
    count[1] = dual_eq_dim_;
    dim_mem[1] = dual_eq_dim_;

    d_mem = H5Screate_simple(2, dim_mem, nullptr);
    hid_t dual_eq = H5Dget_space(ds_dual_eq_);
    H5Sselect_hyperslab(dual_eq, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    H5Dwrite(ds_dual_eq_, H5T_NATIVE_DOUBLE, d_mem, dual_eq, H5P_DATASET_XFER_DEFAULT, state.dual_eq.data());

    H5Sclose(d_mem);
    H5Sclose(dual_eq);
  }

  index_++;
}
