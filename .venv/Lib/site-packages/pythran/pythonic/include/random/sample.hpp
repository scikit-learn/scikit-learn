#ifndef PYTHONIC_INCLUDE_RANDOM_SAMPLE_HPP
#define PYTHONIC_INCLUDE_RANDOM_SAMPLE_HPP

#include "pythonic/include/random/random.hpp"
#include "pythonic/include/utils/functor.hpp"

#include "pythonic/include/types/list.hpp"

PYTHONIC_NS_BEGIN

namespace random
{
  template <class Iterable>
  types::list<typename std::iterator_traits<
      typename std::remove_cv_t<std::remove_reference_t<Iterable>>::iterator>::value_type>
  sample(Iterable &&s, size_t k);

  DEFINE_FUNCTOR(pythonic::random, sample);
} // namespace random
PYTHONIC_NS_END

#endif
