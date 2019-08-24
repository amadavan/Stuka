
workspace(name = "stuka")

# Local configuration
load("@//third_party/systemlibs:syslibs_configure.bzl", "syslibs_configure")
load("@//third_party/py:python_configure.bzl", "python_configure")
load("@//third_party:repo.bzl", "stuka_http_archive")

syslibs_configure(name = "local_config_syslibs")
python_configure(name = "local_config_python")

# Local libraries
load("@//third_party/gurobi:build_defs.bzl", "gurobi_repository")

gurobi_repository(
        name = "gurobi",
        build_file = "@//third_party/gurobi:gurobi.BUILD",
)

# External libraries
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

# Rule repository
http_archive(
   name = "rules_foreign_cc",
   strip_prefix = "rules_foreign_cc-master",
   sha256 = "8ffa95260f5241f501f50ac77ef36863fd610ca1cecb35646d656296f3e4bc82",
   urls = [
       "https://github.com/bazelbuild/rules_foreign_cc/archive/master.zip",
       ]
)

load("@rules_foreign_cc//:workspace_definitions.bzl", "rules_foreign_cc_dependencies")
rules_foreign_cc_dependencies()

# Libraries
stuka_http_archive(
    name = "openblas",
    build_file = "@//third_party:openblas.BUILD",
    system_build_file = "@//third_party/systemlibs:openblas.BUILD",
    sha256 = "e64c8fe083832ffbc1459ab6c72f71d53afd3b36e8497c922a15a06b72e9002f",
    strip_prefix = "OpenBLAS-0.3.6",
    urls = [
        "https://github.com/xianyi/OpenBLAS/archive/v0.3.6.tar.gz",
    ],
)

http_archive(
    name = "blis",
    build_file = "@//third_party:blis.BUILD",
    sha256 = "ad5765cc3f492d0c663f494850dafc4d72f901c332eb442f404814ff2995e5a9",
    strip_prefix = "blis-0.6.0",
    urls = [
        "https://github.com/flame/blis/archive/0.6.0.tar.gz"
    ],
)

http_archive(
    name = "flame",
    build_file = "@//third_party:flame.BUILD",
    sha256 = "e7189b750890bd781fe773f366b374518dd1d89a6513d3d6261bf549826384d1",
    strip_prefix = "libflame-5.1.0",
    urls = [
        "https://github.com/flame/libflame/archive/5.1.0.tar.gz",
    ],
)

http_archive(
    name = "eigen",
    build_file = "@//third_party:eigen.BUILD",
    sha256 = "7e84ef87a07702b54ab3306e77cea474f56a40afa1c0ab245bb11725d006d0da",
    strip_prefix = "eigen-eigen-323c052e1731",
    urls = [
        "http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz",
    ],
)

http_archive(
    name = "suitesparse",
    build_file = "@//third_party:suitesparse.BUILD",
    sha256 = "374dd136696c653e34ef3212dc8ab5b61d9a67a6791d5ec4841efb838e94dbd1",
    strip_prefix = "SuiteSparse",
    urls = [
        "http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.4.0.tar.gz"
    ],
)

stuka_http_archive(
    name = "hdf5",
    build_file = "@//third_party:hdf5.BUILD",
    system_build_file = "@//third_party/systemlibs:hdf5.BUILD",
    sha256 = "6d4ce8bf902a97b050f6f491f4268634e252a63dadd6656a1a9be5b7b7726fa8",
    strip_prefix = "hdf5-1.10.5",
    urls = [
        "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz",
        ]
)

http_archive(
    name = "pybind11",
    build_file = "@//third_party:pybind11.BUILD",
    sha256 = "b69e83658513215b8d1443544d0549b7d231b9f201f6fc787a2b2218b408181e",
    strip_prefix = "pybind11-2.2.4",
    urls = [
        "https://github.com/pybind/pybind11/archive/v2.2.4.tar.gz",
    ],
)

stuka_http_archive(
    name = "zlib",
    build_file = "@//third_party:zlib.BUILD",
    system_build_file = "@//third_party/systemlibs:zlib.BUILD",
    sha256 = "c3e5e9fdd5004dcb542feda5ee4f0ff0744628baf8ed2dd5d66f8ca1197cb1a1",
    strip_prefix = "zlib-1.2.11",
#    system_build_file = "@//third_party/systemlibs:zlib.BUILD",
    urls = [
        "http://mirror.tensorflow.org/zlib.net/zlib-1.2.11.tar.gz",
        "https://zlib.net/zlib-1.2.11.tar.gz",
    ],
)

stuka_http_archive(
    name = "szip",
    build_file = "@//third_party:szip.BUILD",
    system_build_file = "@//third_party/systemlibs:szip.BUILD",
    sha256 = "21ee958b4f2d4be2c9cabfa5e1a94877043609ce86fde5f286f105f7ff84d412",
    strip_prefix = "szip-2.1.1",
    urls = [
        "https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz",
    ]
)

load("@//third_party/mkl:build_defs.bzl", "mkl_repository")

mkl_repository(
    name = "mkl_linux",
    build_file = "@//third_party/mkl:mkl.BUILD",
    sha256 = "f4129843d5c2996419f96f10928edd02b2150998861a088dc7cfa1b6a058102a",
    strip_prefix = "mklml_lnx_2019.0.3.20190220",
    urls = [
        "http://mirror.tensorflow.org/github.com/intel/mkl-dnn/releases/download/v0.18/mklml_lnx_2019.0.3.20190220.tgz",
        "https://github.com/intel/mkl-dnn/releases/download/v0.18/mklml_lnx_2019.0.3.20190220.tgz",
    ],
)

mkl_repository(
    name = "mkl_windows",
    build_file = "@//third_party/mkl:mkl.BUILD",
    sha256 = "eae0c49a7ed738f0ed97b897e952eaa881feddfa665017a8d5d9d79fd38964b4",
    strip_prefix = "mklml_win_2019.0.3.20190220",
    urls = [
        "http://mirror.tensorflow.org/github.com/intel/mkl-dnn/releases/download/v0.18/mklml_win_2019.0.3.20190220.zip",
        "https://github.com/intel/mkl-dnn/releases/download/v0.18/mklml_win_2019.0.3.20190220.zip",
    ],
)

mkl_repository(
    name = "mkl_darwin",
    build_file = "@//third_party/mkl:mkl.BUILD",
    sha256 = "53fdcd7e31c309bb6af869d82987d9c6414c1b957d63d10a9caa9ad077643d99",
    strip_prefix = "mklml_mac_2019.0.3.20190220",
    urls = [
        "http://mirror.tensorflow.org/github.com/intel/mkl-dnn/releases/download/v0.18/mklml_mac_2019.0.3.20190220.tgz",
        "https://github.com/intel/mkl-dnn/releases/download/v0.18/mklml_mac_2019.0.3.20190220.tgz",
    ],
)