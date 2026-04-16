#ifndef PYTHONIC_NUMPY_COMPLEX_HPP
#define PYTHONIC_NUMPY_COMPLEX_HPP

#include "pythonic/include/numpy/complex.hpp"

#include "pythonic/types/complex.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {
    inline std::complex<double> complex(double v, double v2)
    {
      return {v, v2};
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME complex
#define NUMPY_NARY_FUNC_SYM details::complex
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
