//
// Created by Avinash Madavan on 2019-05-08.
//

#ifndef STUKA_CALLBACK_SAVE_HDF5_H
#define STUKA_CALLBACK_SAVE_HDF5_H

#include <hdf5.h>
//#include <H5Cpp.h>

#include "base_callback.h"

#define HDF5_CHUNK_CACHE 1000 * 125 * 4

namespace stuka { namespace util { namespace callback {
  class SaveHDF5 : public BaseCallback {
  public:
    SaveHDF5(const std::string &filename, const size_t n_entries = 0, const bool compress = false);

    ~SaveHDF5() final;

    void callback(const OptimizeState state) override;

    void initialize(const OptimizeState state) override;

    void finish(const OptimizeState state) override;

  private:
    hid_t file_, ds_x_, ds_dual_ub_, ds_dual_eq_, ds_nit_;
//    H5::H5File file_;
//    H5::DataSet ds_x_, ds_dual_ub_, ds_dual_eq_, ds_nit_;
    size_t n_p_, n_entries_, index_;

    size_t x_dim_, dual_ub_dim_, dual_eq_dim_;

    bool compress_;

    void save(const OptimizeState state);
  };
}}}

#endif //STUKA_CALLBACK_SAVE_HDF5_H