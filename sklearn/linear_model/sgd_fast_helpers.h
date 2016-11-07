// We cannot directly reuse the npy_isfinite from npy_math.h as numpy
// and scikit-learn are not necessarily built with the same compiler.
#ifdef _MSC_VER
# include <float.h>
# define skl_isfinite _finite
#else
# include <numpy/npy_math.h>
# define skl_isfinite npy_isfinite
#endif
