#ifndef PYTHONIC_SCIPY_SPECIAL_YV_HPP
#define PYTHONIC_SCIPY_SPECIAL_YV_HPP

#include "pythonic/include/scipy/special/yv.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

#include "pythonic/utils/boost_local_config.hpp"
#include <boost/math/special_functions/bessel.hpp>

PYTHONIC_NS_BEGIN

namespace scipy
{
  namespace special
  {
    namespace details
    {
      template <class T0, class T1>
      double yv(T0 x, T1 y)
      {
        using namespace boost::math::policies;
        return boost::math::cyl_neumann(x, y, make_policy(promote_double<true>()));
      }
    } // namespace details

#define NUMPY_NARY_FUNC_NAME yv
#define NUMPY_NARY_FUNC_SYM details::yv
#include "pythonic/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
