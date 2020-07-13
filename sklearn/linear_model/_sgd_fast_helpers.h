// We cannot directly reuse the npy_isfinite from npy_math.h as numpy
// and scikit-learn are not necessarily built with the same compiler.
// When re-declaring the functions in the template for cython
// specific for each parameter input type, it needs to be 2 different functions
// as cython doesn't support function overloading.
#ifdef _MSC_VER
# include <float.h>
# define skl_isfinite _finite
# define skl_isfinite32 _finite
# define skl_isfinite64 _finite
#else
# include <numpy/npy_math.h>
# define skl_isfinite npy_isfinite
# define skl_isfinite32 npy_isfinite
# define skl_isfinite64 npy_isfinite
#endif
