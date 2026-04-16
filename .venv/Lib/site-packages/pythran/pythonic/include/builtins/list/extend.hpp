#ifndef PYTHONIC_INCLUDE_BUILTIN_LIST_EXTEND_HPP
#define PYTHONIC_INCLUDE_BUILTIN_LIST_EXTEND_HPP

#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace list
  {

    template <class T0, class T1>
    std::enable_if_t<!std::is_same<std::decay_t<T0>, types::empty_list>::value, types::none_type>
    extend(T0 &&seq, T1 const &add);

    template <class T0, class T1>
    std::enable_if_t<std::is_same<std::decay_t<T0>, types::empty_list>::value, types::none_type>
    extend(T0 &&seq, T1 const &add);

    DEFINE_FUNCTOR(pythonic::builtins::list, extend);
  } // namespace list
} // namespace builtins
PYTHONIC_NS_END
#endif
