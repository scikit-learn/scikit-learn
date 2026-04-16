#ifndef PYTHONIC_INCLUDE_ITERTOOLS_IFILTER_HPP
#define PYTHONIC_INCLUDE_ITERTOOLS_IFILTER_HPP

#include "pythonic/include/builtins/filter.hpp"
#include "pythonic/include/itertools/common.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/iterator.hpp"

#include <iterator>
#include <type_traits>

PYTHONIC_NS_BEGIN

namespace itertools
{

  template <typename Operator, typename List0>
  details::filter<std::remove_cv_t<std::remove_reference_t<Operator>>,
                  std::remove_cv_t<std::remove_reference_t<List0>>>
  ifilter(Operator &&_op, List0 &&_seq);

  DEFINE_FUNCTOR(pythonic::itertools, ifilter);
} // namespace itertools

PYTHONIC_NS_END

#endif
