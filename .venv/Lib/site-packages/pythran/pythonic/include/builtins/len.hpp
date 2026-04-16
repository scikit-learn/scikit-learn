#ifndef PYTHONIC_INCLUDE_BUILTIN_LEN_HPP
#define PYTHONIC_INCLUDE_BUILTIN_LEN_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/yield.hpp"

#include <iterator>
#include <tuple>

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class... Types>
  long len(std::tuple<Types...> const &);

  template <class T>
  std::enable_if_t<types::has_size<T>::value, long> len(T const &t);

  DEFINE_FUNCTOR(pythonic::builtins, len);
} // namespace builtins
PYTHONIC_NS_END
#endif
