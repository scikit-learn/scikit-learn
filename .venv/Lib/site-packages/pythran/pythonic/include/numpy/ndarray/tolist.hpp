#ifndef PYTHONIC_INCLUDE_NUMPY_NDARRAY_TOLIST_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDARRAY_TOLIST_HPP

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace ndarray
  {
    template <class T, size_t N>
    struct tolist_type {
      using type = types::list<typename tolist_type<T, N - 1>::type>;
    };

    template <class T>
    struct tolist_type<T, 1> {
      using type = types::list<T>;
    };

    template <class T, class pS>
    std::enable_if_t<std::tuple_size<pS>::value == 1, types::list<T>>
    tolist(types::ndarray<T, pS> const &expr);

    template <class T, class pS>
    std::enable_if_t<std::tuple_size<pS>::value != 1,
                     typename tolist_type<T, std::tuple_size<pS>::value>::type>
    tolist(types::ndarray<T, pS> const &expr);

    NUMPY_EXPR_TO_NDARRAY0_DECL(tolist);
    DEFINE_FUNCTOR(pythonic::numpy::ndarray, tolist);
  } // namespace ndarray
} // namespace numpy
PYTHONIC_NS_END

#endif
