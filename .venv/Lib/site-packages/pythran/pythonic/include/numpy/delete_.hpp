#ifndef PYTHONIC_INCLUDE_NUMPY_DELETE_HPP
#define PYTHONIC_INCLUDE_NUMPY_DELETE_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, types::pshape<long>> delete_(types::ndarray<T, pS> const &a, long index,
                                                 types::none_type axis = builtins::None);

  template <class T, class pS, class I>
  std::enable_if_t<!std::is_scalar<I>::value, types::ndarray<T, types::pshape<long>>>
  delete_(types::ndarray<T, pS> const &in, I const &indices,
          types::none_type axis = builtins::None);

  NUMPY_EXPR_TO_NDARRAY0_DECL(delete_);
  DEFINE_FUNCTOR(pythonic::numpy, delete_);
} // namespace numpy
PYTHONIC_NS_END

#endif
