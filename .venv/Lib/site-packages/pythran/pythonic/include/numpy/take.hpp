#ifndef PYTHONIC_INCLUDE_NUMPY_TAKE_HPP
#define PYTHONIC_INCLUDE_NUMPY_TAKE_HPP

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class F, class T>
  auto take(T &&expr, F &&indices) -> decltype(std::forward<T>(expr)[std::forward<F>(indices)]);

  DEFINE_FUNCTOR(pythonic::numpy, take);
} // namespace numpy
PYTHONIC_NS_END

#endif
