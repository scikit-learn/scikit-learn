#ifndef PYTHONIC_INCLUDE_BUILTIN_ISINSTANCE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_ISINSTANCE_HPP

#include "pythonic/builtins/pythran/is_none.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"

PYTHONIC_NS_BEGIN
namespace types
{
  class str;

  template <class Ty0, class Ty1>
  struct isinstance : std::conditional<std::is_same<Ty0, Ty1>::value, true_type, false_type> {
  };

  // some specialization

  template <>
  struct isinstance<char const *, str> {
    using type = true_type;
  };
  template <>
  struct isinstance<str, char const *> {
    using type = true_type;
  };
} // namespace types

namespace builtins
{
  namespace details
  {
    template <class Obj, class Cls>
    struct isinstance {
      using type =
          typename types::isinstance<Obj, std::decay_t<decltype(std::declval<Cls>()())>>::type;
    };

    template <class Obj, class... Clss>
    struct isinstance<Obj, std::tuple<Clss...>> {
      using type = std::conditional_t<
          utils::any_of<std::is_same<
              typename types::isinstance<Obj, std::decay_t<decltype(std::declval<Clss>()())>>::type,
              types::true_type>::value...>::value,
          types::true_type, types::false_type>;
    };
  } // namespace details

  template <class Obj, class Cls>
  typename details::isinstance<Obj, Cls>::type isinstance(Obj, Cls)
  {
    return {};
  }

  DEFINE_FUNCTOR(pythonic::builtins, isinstance);
} // namespace builtins
PYTHONIC_NS_END

#endif
