#ifndef PYTHONIC_NUMPY_INSERT_HPP
#define PYTHONIC_NUMPY_INSERT_HPP

#include "pythonic/include/numpy/insert.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/traits.hpp"
#include "pythonic/utils/functor.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class T, class pS, class I, class F>
  std::enable_if_t<types::is_iterable<I>::value && types::is_iterable<F>::value,
                   types::ndarray<T, types::pshape<long>>>
  insert(types::ndarray<T, pS> in, I const &indices, F const &data, types::none_type axis)
  {
    types::ndarray<T, types::pshape<long>> out(
        types::pshape<long>(long(in.flat_size() + std::min(indices.flat_size(), data.flat_size()))),
        builtins::None);
    auto out_iter = out.fbegin();
    auto in_iter = in.fbegin();
    auto data_iter = data.begin();
    for (long index : indices) {
      out_iter = std::copy(in_iter, in.fbegin() + index, out_iter);
      *out_iter++ = *data_iter++;
      in_iter = in.fbegin() + index;
    }
    std::copy(in_iter, in.fend(), out_iter);
    return out;
  }

  template <class T, class pS, class I, class F>
  std::enable_if_t<types::is_iterable<I>::value && !types::is_iterable<F>::value,
                   types::ndarray<T, types::pshape<long>>>
  insert(types::ndarray<T, pS> in, I const &indices, F const &data, types::none_type axis)
  {
    return insert(in, indices, types::list<F>({data}), axis);
  }

  template <class T, class pS, class I, class F>
  std::enable_if_t<!types::is_iterable<I>::value && types::is_iterable<F>::value,
                   types::ndarray<T, types::pshape<long>>>
  insert(types::ndarray<T, pS> in, I const &indices, F const &data, types::none_type axis)
  {
    return insert(in, types::list<I>({indices}), {data}, axis);
  }

  template <class T, class pS, class I, class F>
  std::enable_if_t<!types::is_iterable<I>::value && !types::is_iterable<F>::value,
                   types::ndarray<T, types::pshape<long>>>
  insert(types::ndarray<T, pS> in, I const &indices, F const &data, types::none_type axis)
  {
    return insert(in, types::list<I>({indices}), types::list<F>({data}), axis);
  }

  template <class E, class... Args>
  E insert(E, Args const &...)
  {
    throw std::runtime_error("insert only partially supported");
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
