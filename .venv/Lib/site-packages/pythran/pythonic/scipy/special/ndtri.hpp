#ifndef PYTHONIC_SCIPY_SPECIAL_NDTRI_HPP
#define PYTHONIC_SCIPY_SPECIAL_NDTRI_HPP

#include "pythonic/include/scipy/special/ndtri.hpp"

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
      double ndtri(T x)
      {
        using namespace boost::math::policies;
        boost::math::normal dist(0.0, 1.0);
        return quantile(dist, x);
      }
    } // namespace details

#define NUMPY_NARY_FUNC_NAME ndtri
#define NUMPY_NARY_FUNC_SYM details::ndtri
#include "pythonic/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
