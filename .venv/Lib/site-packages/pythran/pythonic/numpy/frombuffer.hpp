#ifndef PYTHONIC_NUMPY_FROMBUFFER_HPP
#define PYTHONIC_NUMPY_FROMBUFFER_HPP

#include "pythonic/include/numpy/frombuffer.hpp"

#include "pythonic/types/list.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

#include <limits>
#include <sstream>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype>
  types::ndarray<typename dtype::type, types::pshape<long>>
  frombuffer(types::str const &string, dtype d, long count, long offset)
  {
    if (count < 0)
      count = string.size() / sizeof(typename dtype::type);
    types::pshape<long> shape = count;
    utils::shared_ref<types::raw_array<typename dtype::type>> buffer(std::get<0>(shape));
    auto const *tstring = reinterpret_cast<typename dtype::type const *>(string.c_str()) + offset;
    std::copy(tstring, tstring + std::get<0>(shape), buffer->data);
    return {buffer, shape};
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
