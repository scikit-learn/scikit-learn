#ifndef PYTHONIC_INCLUDE_RANDOM_SHUFFLE_HPP
#define PYTHONIC_INCLUDE_RANDOM_SHUFFLE_HPP

#include "pythonic/include/random/random.hpp"
#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace random
{
  template <class T>
  types::none_type shuffle(T &seq);
  template <class T, class function>
  types::none_type shuffle(T &seq, function &&randf);

  DEFINE_FUNCTOR(pythonic::random, shuffle)
} // namespace random
PYTHONIC_NS_END

#endif
