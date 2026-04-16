#ifndef PYTHONIC_INCLUDE_NUMPY_COMPLEX_HPP
#define PYTHONIC_INCLUDE_NUMPY_COMPLEX_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {

    std::complex<double> complex(double v = 0, double v2 = 0.);
  }

#define NUMPY_NARY_FUNC_NAME complex
#define NUMPY_NARY_FUNC_SYM details::complex
#define NUMPY_NARY_EXTRA_METHOD using type = std::complex<double>;
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
