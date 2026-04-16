#ifndef PYTHONIC_NUMPY_FROMITER_HPP
#define PYTHONIC_NUMPY_FROMITER_HPP

#include "pythonic/include/numpy/fromiter.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class Iterable, class dtype>
  types::ndarray<typename std::remove_cv_t<std::remove_reference_t<Iterable>>::value_type,
                 types::pshape<long>>
  fromiter(Iterable &&iterable, dtype d, long count)
  {
    using T = typename std::remove_cv_t<std::remove_reference_t<Iterable>>::value_type;
    if (count < 0) {
      types::list<T> buffer(0);
      std::copy(iterable.begin(), iterable.end(), std::back_inserter(buffer));
      return {buffer};
    } else {
      utils::shared_ref<types::raw_array<T>> buffer(count);
      std::copy_n(iterable.begin(), count, buffer->data);
      types::array_tuple<long, 1> shape = {count};
      return {buffer, shape};
    }
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
