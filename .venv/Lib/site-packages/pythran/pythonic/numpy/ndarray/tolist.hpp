#ifndef PYTHONIC_NUMPY_NDARRAY_TOLIST_HPP
#define PYTHONIC_NUMPY_NDARRAY_TOLIST_HPP

#include "pythonic/include/numpy/ndarray/tolist.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace ndarray
  {

    template <class T, class pS>
    std::enable_if_t<std::tuple_size<pS>::value == 1, types::list<T>>
    tolist(types::ndarray<T, pS> const &expr)
    {
      return {expr.fbegin(), expr.fend()};
    }

    template <class T, class pS>
    std::enable_if_t<std::tuple_size<pS>::value != 1,
                     typename tolist_type<T, std::tuple_size<pS>::value>::type>
    tolist(types::ndarray<T, pS> const &expr)
    {
      typename tolist_type<T, std::tuple_size<pS>::value>::type out(0);
      for (auto const &elts : expr)
        out.push_back(tolist(elts));
      return out;
    }

    NUMPY_EXPR_TO_NDARRAY0_IMPL(tolist);
  } // namespace ndarray
} // namespace numpy
PYTHONIC_NS_END

#endif
