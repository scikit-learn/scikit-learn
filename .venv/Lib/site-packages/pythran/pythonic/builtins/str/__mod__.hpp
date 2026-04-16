#ifndef PYTHONIC_BUILTIN_STR_MOD_HPP
#define PYTHONIC_BUILTIN_STR_MOD_HPP

#include "pythonic/builtins/str/__mod__.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

#include <boost/format.hpp>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    namespace details
    {
      template <class Tuple>
      void fmt(boost::format &f, Tuple const &a, utils::int_<1>)
      {
        f % std::get<std::tuple_size<Tuple>::value - 1>(a);
      }
      template <class Tuple, size_t I>
      void fmt(boost::format &f, Tuple const &a, utils::int_<I>)
      {
        fmt(f % std::get<std::tuple_size<Tuple>::value - I>(a), a, utils::int_<I - 1>());
      }
    } // namespace details

    template <class T>
    types::str __mod__(types::str const &s, T const &arg)
    {
      const boost::format fmter(s.chars());
      return (boost::format(fmter) % arg).str();
    }

    template <class... Ts>
    types::str __mod__(types::str const &s, std::tuple<Ts...> const &args)
    {
      boost::format fmter(s.chars());
      details::fmt(fmter, args, utils::int_<sizeof...(Ts)>());
      return fmter.str();
    }
    template <size_t N, class T>
    types::str __mod__(types::str const &s, types::array_tuple<T, N> const &args)
    {
      boost::format fmter(s.chars());
      details::fmt(fmter, args, utils::int_<N>());
      return fmter.str();
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
