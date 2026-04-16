#ifndef PYTHONIC_RANDOM_SAMPLE_HPP
#define PYTHONIC_RANDOM_SAMPLE_HPP

#include "pythonic/include/random/sample.hpp"

#include "pythonic/random/random.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/functor.hpp"

#include "pythonic/types/list.hpp"

PYTHONIC_NS_BEGIN

namespace random
{
  template <class Iterable>
  types::list<typename std::iterator_traits<
      typename std::remove_cv_t<std::remove_reference_t<Iterable>>::iterator>::value_type>
  sample(Iterable &&s, size_t k)
  {
    using value_type = typename std::iterator_traits<
        typename std::remove_cv_t<std::remove_reference_t<Iterable>>::iterator>::value_type;
    types::list<value_type> tmp(s.begin(), s.end());
    std::vector<size_t, utils::allocator<size_t>> indices(tmp.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::random_shuffle(indices.begin(), indices.end());
    types::list<value_type> out(k);
    for (size_t i = 0; i < k; i++)
      out[i] = tmp[indices[i]];
    return out;
  }
} // namespace random
PYTHONIC_NS_END

#endif
