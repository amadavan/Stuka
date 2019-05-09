//
// Created by Avinash Madavan on 2019-05-08.
//

#ifndef STUKA_CALLBACK_SAVE_HDF5_H
#define STUKA_CALLBACK_SAVE_HDF5_H

#include "H5Cpp.h"

#include "base_callback.h"

namespace stuka { namespace util { namespace callback {
  class SaveHDF5 : public BaseCallback {
  public:
    SaveHDF5(const std::string &filename, const size_t n_entries = 0);

    void callback(const OptimizeState state) override;

    void initialize(const OptimizeState state) override;

  private:
    H5::H5File file_;
    H5::DataSet ds_x_, ds_dual_ub_, ds_dual_eq_;
    size_t n_entries_;

    hsize_t x_dim_, dual_ub_dim_, dual_eq_dim_;
  };
}}}

#endif //STUKA_CALLBACK_SAVE_HDF5_H