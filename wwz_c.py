import numpy as np
import ctypes
import numpy.ctypeslib

# below should be modified
_lib_wwz = np.ctypeslib.load_library("libwwz.so", "/Users/kawana/anaconda3/envs/myenv/lib/python3.6/site-packages/wwz")

#  _DOUBLE_PP = np.ctypeslib.ndpointer(dtype = np.uintp, ndim=1, flags="C")
_ndtype_1darray = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")

_lib_wwz.wwz_c.argtypes = [
        _ndtype_1darray,
        _ndtype_1darray,
        _ndtype_1darray,
        _ndtype_1darray,
        ctypes.c_double,
        _ndtype_1darray,
        _ndtype_1darray,
        #  _DOUBLE_PP,
        #  _DOUBLE_PP,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int
        ]

_lib_wwz.wwz_c.restype = None


def wwz(t_sample_arr,
        freq_sample_arr,
        t_data_arr,
        x_data_arr,
        c_coef = 0.0125):
    """
    return WWZ(tau, freq) (kind of power) and WWA(tau, freq) (amplitude)

    -----Arguments-----
    ---- Argument -- ------- dtype ------- -----values --------------
    t_sample_arr     numpy ndarray(ndim=1) array of time sampled
    freq_sample_arr  numpy ndarray(ndim=1) array of frequency sampled
    t_data_arr       numpy ndarray(ndim=1) array of time of data
    x_data_arr       numpy ndarray(ndim=1) array of value of data
    c_coef(optional) float                 coefficeint. dafault value is 0.0125

    -----Return values-----
    (wwz, wwa)
    wwz: np.ndarray(dtype=np.float64, ndim=2) shape == (t_sample_arr, freq_sample_arr)
    wwa: np.ndarray(dtype=np.float64, ndim=2) shape == (t_sample_arr, freq_sample_arr)

    """

    omega_sample_arr = 2 * np.pi * freq_sample_arr

    wwa_arr = np.zeros((t_sample_arr.size * omega_sample_arr.size), dtype=np.float64)
    wwz_arr = np.zeros((t_sample_arr.size * omega_sample_arr.size), dtype=np.float64)

    _lib_wwz.wwz_c(
            t_sample_arr.astype(np.float64, copy=False),
            omega_sample_arr.astype(np.float64, copy=False),
            t_data_arr.astype(np.float64, copy=False),
            x_data_arr.astype(np.float64, copy=False),
            c_coef,
            wwz_arr,
            wwa_arr,
            t_sample_arr.size,
            omega_sample_arr.size,
            t_data_arr.size
            )

    return wwz_arr.reshape(t_sample_arr.size, omega_sample_arr.size), wwa_arr.reshape(t_sample_arr.size, omega_sample_arr.size)

#  cdef extern from "wwz_my.c":
#      cdef void wwz_c(double t_sample_arr[],
#                      double omega_sample_arr[],
#                      double t_data_arr[],
#                      double x_data_arr[],
#                      double c_coef,
#                      double **out_wwz,
#                      double **out_wwa,
#                      int    N_t_sample_arr,
#                      int    N_omega_sample_arr,
#                      int    N_t_data_arr)
#
#  DTYPE = np.float64
#  ctypedef np.float64_t DTYPE_t


