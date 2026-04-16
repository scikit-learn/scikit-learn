#ifndef PYTHONIC_SCIPY_SPECIAL_HANKEL1_HPP
#define PYTHONIC_SCIPY_SPECIAL_HANKEL1_HPP

#include "pythonic/include/scipy/special/hankel1.hpp"

#include "pythonic/types/complex.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

#include "pythonic/utils/boost_local_config.hpp"
#include <boost/math/special_functions/hankel.hpp>

PYTHONIC_NS_BEGIN

namespace scipy
{
  namespace special
  {
    namespace details
    {
      template <class T0, class T1>
      std::complex<double> hankel1(T0 x, T1 y)
      {
        return boost::math::cyl_hankel_1(x, y);
      }
    } // namespace details

#define NUMPY_NARY_FUNC_NAME hankel1
#define NUMPY_NARY_FUNC_SYM details::hankel1
#include "pythonic/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
