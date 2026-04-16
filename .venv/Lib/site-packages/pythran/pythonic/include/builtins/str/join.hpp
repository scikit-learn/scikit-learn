#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_JOIN_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_JOIN_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    /* Join for string.join(string) */
    template <class S>
    types::str join(S const &s, types::str const &iterable);

    /* Join for string.join(random acces iter but ! on string) */
    template <class S, class Iterable>
    std::enable_if_t<
        !std::is_same<std::remove_cv_t<std::remove_reference_t<Iterable>>, types::str>::value &&
            std::is_same<typename std::iterator_traits<typename std::remove_reference_t<
                             Iterable>::iterator>::iterator_category,
                         std::random_access_iterator_tag>::value,
        types::str>
    join(S const &s, Iterable &&iterable);

    /* Join for string.join(forward iterator) */
    template <class S, class Iterable>
    std::enable_if_t<!std::is_same<typename std::iterator_traits<typename std::remove_reference_t<
                                       Iterable>::iterator>::iterator_category,
                                   std::random_access_iterator_tag>::value,
                     types::str>
    join(S const &s, Iterable &&iterable);

    DEFINE_FUNCTOR(pythonic::builtins::str, join);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
