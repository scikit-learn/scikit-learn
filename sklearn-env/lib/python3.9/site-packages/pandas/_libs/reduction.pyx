import numpy as np

cimport numpy as cnp

cnp.import_array()

from pandas._libs.util cimport is_array


cdef cnp.dtype _dtype_obj = np.dtype("object")


cpdef check_result_array(object obj, object dtype):
    # Our operation is supposed to be an aggregation/reduction. If
    #  it returns an ndarray, this likely means an invalid operation has
    #  been passed. See test_apply_without_aggregation, test_agg_must_agg
    if is_array(obj):
        if dtype != _dtype_obj:
            # If it is object dtype, the function can be a reduction/aggregation
            #  and still return an ndarray e.g. test_agg_over_numpy_arrays
            raise ValueError("Must produce aggregated value")


cpdef inline extract_result(object res):
    """ extract the result object, it might be a 0-dim ndarray
        or a len-1 0-dim, or a scalar """
    if hasattr(res, "_values"):
        # Preserve EA
        res = res._values
        if res.ndim == 1 and len(res) == 1:
            # see test_agg_lambda_with_timezone, test_resampler_grouper.py::test_apply
            res = res[0]
    return res
