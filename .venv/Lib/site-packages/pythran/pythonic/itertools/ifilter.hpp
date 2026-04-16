#ifndef PYTHONIC_ITERTOOLS_IFILTER_HPP
#define PYTHONIC_ITERTOOLS_IFILTER_HPP

#include "pythonic/builtins/filter.hpp"
#include "pythonic/include/itertools/ifilter.hpp"

PYTHONIC_NS_BEGIN

namespace itertools
{

  template <typename Operator, typename List0>
  details::filter<std::remove_cv_t<std::remove_reference_t<Operator>>,
                  std::remove_cv_t<std::remove_reference_t<List0>>>
  ifilter(Operator &&_op, List0 &&_seq)
  {
    return {std::forward<Operator>(_op), std::forward<List0>(_seq)};
  }
} // namespace itertools
PYTHONIC_NS_END
#endif
