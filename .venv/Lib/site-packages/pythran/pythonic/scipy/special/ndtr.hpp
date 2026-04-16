#ifndef PYTHONIC_SCIPY_SPECIAL_NDTR_HPP
#define PYTHONIC_SCIPY_SPECIAL_NDTR_HPP

#include "pythonic/include/scipy/special/ndtr.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

#include "pythonic/utils/boost_local_config.hpp"
#include <boost/math/distributions/normal.hpp>

PYTHONIC_NS_BEGIN

namespace scipy
{
  namespace special
  {
    namespace details
    {
      template <class T>
      double ndtr(T x)
      {
        using namespace boost::math::policies;
        boost::math::normal dist(0.0, 1.0);
        return cdf(dist, x);
      }
    } // namespace details

#define NUMPY_NARY_FUNC_NAME ndtr
#define NUMPY_NARY_FUNC_SYM details::ndtr
#include "pythonic/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
