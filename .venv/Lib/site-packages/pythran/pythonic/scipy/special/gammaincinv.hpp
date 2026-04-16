#ifndef PYTHONIC_SCIPY_SPECIAL_GAMMAINCINV_HPP
#define PYTHONIC_SCIPY_SPECIAL_GAMMAINCINV_HPP

#include "pythonic/include/scipy/special/gammaincinv.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

#include "pythonic/utils/boost_local_config.hpp"
#include <boost/math/special_functions/gamma.hpp>

PYTHONIC_NS_BEGIN

namespace scipy
{
  namespace special
  {
    namespace details
    {
      template <class T0, class T1>
      double gammaincinv(T0 a, T1 p)
      {
        using namespace boost::math::policies;
        return boost::math::gamma_p_inv(a, p, make_policy(promote_double<true>()));
      }
    } // namespace details

#define NUMPY_NARY_FUNC_NAME gammaincinv
#define NUMPY_NARY_FUNC_SYM details::gammaincinv
#include "pythonic/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
