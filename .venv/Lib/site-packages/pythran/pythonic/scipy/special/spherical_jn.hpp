#ifndef PYTHONIC_SCIPY_SPECIAL_SPHERICAL_JN_HPP
#define PYTHONIC_SCIPY_SPECIAL_SPHERICAL_JN_HPP

#include "pythonic/include/scipy/special/spherical_jn.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

#include "pythonic/utils/boost_local_config.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>

PYTHONIC_NS_BEGIN

namespace scipy
{
  namespace special
  {
    namespace details
    {
      template <class T0, class T1>
      double spherical_jn(T0 v, T1 x, bool derivative)
      {
        assert(v == (long)v && "only supported for integral value as first arg");
        using namespace boost::math::policies;
        if (derivative) {
          return boost::math::sph_bessel_prime(v, x, make_policy(promote_double<true>()));
        } else {
          return boost::math::sph_bessel(v, x, make_policy(promote_double<true>()));
        }
      }
    } // namespace details

#define NUMPY_NARY_FUNC_NAME spherical_jn
#define NUMPY_NARY_FUNC_SYM details::spherical_jn
#include "pythonic/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
