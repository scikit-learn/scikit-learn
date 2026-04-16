#ifndef PYTHONIC_NUMPY_CONCATENATE_HPP
#define PYTHONIC_NUMPY_CONCATENATE_HPP

#include "pythonic/include/numpy/concatenate.hpp"

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/builtins/sum.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    template <size_t N>
    struct concatenate_helper {
      // list version
      template <class Out, class A>
      void operator()(Out &&out, A const &from, long axis) const
      {
        if (axis == 0) {
          auto out_iter = out.begin();
          for (auto &&ifrom : from)
            out_iter = std::copy(ifrom.begin(), ifrom.end(), out_iter);
        } else {
          using iterator_on_from_value = typename A::value_type::const_iterator;
          std::vector<iterator_on_from_value, utils::allocator<iterator_on_from_value>> ifroms;
          for (auto &ifrom : from)
            ifroms.emplace_back(ifrom.begin());

          using iterator_value_type =
              typename std::iterator_traits<iterator_on_from_value>::value_type;
          std::vector<iterator_value_type, utils::allocator<iterator_value_type>> difroms;

          for (auto &&iout : out) {
            difroms.clear();
            for (auto &ifrom : ifroms)
              difroms.emplace_back(*ifrom);
            concatenate_helper<N - 1>()(iout, difroms, axis - 1);
            for (auto &ifrom : ifroms)
              ++ifrom;
          }
        }
      }
      // array version
      template <class Out, class A, size_t... I>
      void operator()(Out &&out, A const &from, long axis, std::index_sequence<I...>) const
      {
        if (axis == 0) {
          auto out_iter = out.begin();
          (void)std::initializer_list<int>{
              (out_iter = std::copy(std::get<I>(from).begin(), std::get<I>(from).end(), out_iter),
               1)...};
        } else {
          types::array_tuple<typename A::value_type::const_iterator, sizeof...(I)> ifroms = {
              std::get<I>(from).begin()...};

          for (auto &&iout : out) {
            types::array_tuple<
                typename std::iterator_traits<typename A::value_type::const_iterator>::value_type,
                sizeof...(I)>
                difroms = {*std::get<I>(ifroms)...};
            concatenate_helper<N - 1>()(iout, difroms, axis - 1, std::index_sequence<I...>{});
            (void)std::initializer_list<int>{(++std::get<I>(ifroms), 0)...};
          }
        }
      }
      // tuple version
      template <class Out, class... Ts, size_t... I>
      void operator()(Out &&out, std::tuple<Ts...> const &from, long axis,
                      std::index_sequence<I...>) const
      {
        if (axis == 0) {
          auto out_iter = out.begin();
          (void)std::initializer_list<int>{
              (out_iter = std::copy(std::get<I>(from).begin(), std::get<I>(from).end(), out_iter),
               1)...};
        } else {
          auto ifroms = std::make_tuple(std::get<I>(from).begin()...);

          for (auto &&iout : out) {
            auto difroms = std::make_tuple(*std::get<I>(ifroms)...);
            concatenate_helper<N - 1>()(iout, difroms, axis - 1, std::index_sequence<I...>{});
            (void)std::initializer_list<int>{(++std::get<I>(ifroms), 0)...};
          }
        }
      }
    };

    template <>
    struct concatenate_helper<0> {
      // list version - sentinel
      template <class Out, class A>
      void operator()(Out &&buffer, A const &from, long axis) const
      {
      }
      // array version
      template <class Out, class E, size_t... I>
      void operator()(Out &&, E const &, long, std::index_sequence<I...>) const
      {
      }
      // tuple version - sentinel
      template <class Out, class... Ts, size_t... I>
      void operator()(Out &&, std::tuple<Ts...> const &, long, std::index_sequence<I...>) const
      {
      }
    };

    template <class A, size_t... I>
    long concatenate_axis_size(A const &from, long axis, std::index_sequence<I...>)
    {
      long sizes[] = {sutils::getshape(std::get<I>(from))[axis]...};
      return std::accumulate(std::begin(sizes), std::end(sizes), 0L, std::plus<long>());
    }
  } // namespace details

  template <class... Types>
  auto concatenate(std::tuple<Types...> const &args, long axis) -> types::ndarray<
      typename __combined<typename std::decay_t<Types>::dtype...>::type,
      types::array_tuple<long, std::tuple_element_t<0, std::tuple<Types...>>::value>>
  {
    auto constexpr N = std::decay_t<decltype(std::get<0>(args))>::value;
    auto shape = sutils::getshape(std::get<0>(args));
    shape[axis] =
        details::concatenate_axis_size(args, axis, std::make_index_sequence<sizeof...(Types)>{});

    types::ndarray<typename __combined<typename std::decay_t<Types>::dtype...>::type,
                   types::array_tuple<long, std::decay_t<decltype(std::get<0>(args))>::value>>
        result{shape, types::none_type{}};
    details::concatenate_helper<N>()(result, args, axis,
                                     std::make_index_sequence<sizeof...(Types)>{});
    return result;
  }

  template <class E, size_t M, class V>
  types::ndarray<typename E::dtype, types::array_tuple<long, E::value>>
  concatenate(types::array_base<E, M, V> const &args, long axis)
  {
    auto constexpr N = E::value;
    auto shape = sutils::getshape(std::get<0>(args));
    shape[axis] = details::concatenate_axis_size(args, axis, std::make_index_sequence<M>{});
    types::ndarray<typename E::dtype, types::array_tuple<long, E::value>> out(shape,
                                                                              types::none_type{});
    details::concatenate_helper<N>()(out, args, axis, std::make_index_sequence<M>{});
    return out;
  }

  template <class E>
  types::ndarray<typename E::dtype, types::array_tuple<long, E::value>>
  concatenate(types::list<E> const &ai, long axis)
  {
    using return_type = types::ndarray<typename E::dtype, types::array_tuple<long, E::value>>;
    auto constexpr N = return_type::value;
    auto shape = sutils::getshape(ai[0]);
    shape[axis] = std::accumulate(ai.begin(), ai.end(), 0L, [axis](long v, E const &from) {
      return v + sutils::getshape(from)[axis];
    });

    return_type out{shape, types::none_type{}};
    details::concatenate_helper<N>()(out, ai, axis);
    return out;
    ;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
