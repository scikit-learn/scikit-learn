#ifndef PYTHONIC_NUMPY_DELETE_HPP
#define PYTHONIC_NUMPY_DELETE_HPP

#include "pythonic/include/numpy/delete_.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, types::pshape<long>> delete_(types::ndarray<T, pS> const &a, long index,
                                                 types::none_type axis)
  {
    types::ndarray<T, types::pshape<long>> out(types::pshape<long>(long(a.flat_size()) - 1),
                                               builtins::None);
    long n = a.flat_size();
    index = std::min(n, index);
    std::copy(a.buffer + index + 1, a.buffer + n,
              std::copy(a.buffer, a.buffer + index, out.buffer));
    return out;
  }

  template <class T, class pS, class I>
  std::enable_if_t<!std::is_scalar<I>::value, types::ndarray<T, types::pshape<long>>>
  delete_(types::ndarray<T, pS> const &in, I const &indices, types::none_type axis)
  {
    types::ndarray<T, types::pshape<long>> out(
        types::pshape<long>(long(in.flat_size()) - indices.flat_size()), builtins::None);
    auto out_iter = out.buffer;
    auto in_iter = in.buffer;
    for (long index : indices) {
      out_iter = std::copy(in_iter, in.buffer + index, out_iter);
      in_iter = in.buffer + index + 1;
    }
    std::copy(in_iter, in.buffer + in.flat_size(), out_iter);
    return out;
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(delete_);
} // namespace numpy
PYTHONIC_NS_END

#endif
