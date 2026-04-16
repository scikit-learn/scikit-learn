#ifndef PYTHONIC_INCLUDE_TYPES_DYNAMIC_TUPLE_HPP
#define PYTHONIC_INCLUDE_TYPES_DYNAMIC_TUPLE_HPP

#include "pythonic/include/types/assignable.hpp"
#include "pythonic/include/types/nditerator.hpp"
#include "pythonic/include/types/traits.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/allocate.hpp"
#include "pythonic/include/utils/int_.hpp"
#include "pythonic/include/utils/nested_container.hpp"
#include "pythonic/include/utils/seq.hpp"
#include "pythonic/include/utils/shared_ref.hpp"

#include <vector>

PYTHONIC_NS_BEGIN

namespace types
{

  template <typename T>
  struct dynamic_tuple {
    using container_type = std::vector<T, utils::allocator<T>>;
    utils::shared_ref<container_type> data;

    using value_type = T;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator = typename container_type::const_iterator;
    using const_iterator = typename container_type::const_iterator;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reverse_iterator = typename container_type::reverse_iterator;
    using const_reverse_iterator = typename container_type::const_reverse_iterator;

    // minimal ndarray interface
    using dtype = typename utils::nested_container_value_type<dynamic_tuple>::type;
    static const size_t value = utils::nested_container_depth<dynamic_tuple>::value;
    static const bool is_vectorizable = true;
    static const bool is_strided = false;

    // flat_size implementation
    template <class E>
    long _flat_size(E const &e, utils::int_<1>) const;
    template <class E, size_t L>
    long _flat_size(E const &e, utils::int_<L>) const;

    long flat_size() const;

    dynamic_tuple() = default;
    dynamic_tuple(dynamic_tuple const &) = default;
    dynamic_tuple(dynamic_tuple &&) = default;
    dynamic_tuple &operator=(dynamic_tuple &&other) = default;
    dynamic_tuple &operator=(dynamic_tuple const &other) = default;

    template <class Iter>
    dynamic_tuple(Iter start, Iter end) : data(start, end)
    {
    }

    dynamic_tuple(std::initializer_list<T> values) : data(values)
    {
    }

    // Iterators.
    const_iterator begin() const noexcept
    {
      return data->begin();
    }

    const_iterator end() const noexcept
    {
      return data->end();
    }

    const_reverse_iterator rbegin() const noexcept
    {
      return data->rbegin();
    }

    const_reverse_iterator rend() const noexcept
    {
      return data->rend();
    }

    // Capacity.
    size_type size() const noexcept
    {
      return data->size();
    }
    constexpr bool empty() const noexcept
    {
      return data->empty();
    }

    intptr_t id() const;

    // Element access.
    const_reference fast(long n) const
    {
      return (*data)[n];
    }
#ifdef USE_XSIMD
    using simd_iterator = const_simd_nditerator<dynamic_tuple>;
    using simd_iterator_nobroadcast = simd_iterator;
    template <class vectorizer>
    simd_iterator vbegin(vectorizer) const;
    template <class vectorizer>
    simd_iterator vend(vectorizer) const;
#endif

    const_reference operator[](size_type __n) const
    {
      return (*data)[__n < 0 ? __n + size() : __n];
    }

    reference operator[](size_type __n)
    {
      return (*data)[__n < 0 ? __n + size() : __n];
    }

    // operator
    bool operator==(dynamic_tuple<T> const &other) const;

    bool operator!=(dynamic_tuple<T> const &other) const;

    bool operator<(dynamic_tuple<T> const &other) const;
    bool operator<=(dynamic_tuple<T> const &other) const;

    bool operator>(dynamic_tuple<T> const &other) const;
    bool operator>=(dynamic_tuple<T> const &other) const;

    dynamic_tuple<T> operator+(dynamic_tuple<T> const &other) const;

    dynamic_tuple operator[](slice const &s) const
    {
      auto ns = s.normalize(size());
      dynamic_tuple res;
      res.data->reserve(ns.size());
      for (auto i = ns.lower, step = ns.step, n = ns.upper; i != n; i += step) {
        res.data->emplace_back(fast(i));
      }
      return res;
    }

    template <long stride>
    dynamic_tuple operator[](cstride_slice<stride> const &s) const
    {
      auto ns = s.normalize(size());
      if (stride == 1)
        return {begin() + ns.lower, begin() + ns.upper};
      else {
        dynamic_tuple res;
        res.data->reserve(ns.size());
        for (auto i = ns.lower, step = ns.step, n = ns.upper; i != n; i += step) {
          res.data->emplace_back(fast(i));
        }
        return res;
      }
    }

    dynamic_tuple operator[](fast_contiguous_slice const &s) const
    {
      auto ns = s.normalize(size());
      return {begin() + ns.lower, begin() + ns.upper};
    }

    using shape_t = typename shape_builder<dynamic_tuple, value>::type;
    template <size_t I>
    auto shape() const -> decltype(details::extract_shape(*this, utils::int_<I>{}))
    {
      return details::extract_shape(*this, utils::int_<I>{});
    }

    template <class E, size_t N, class S>
    operator array_base<E, N, S>() const
    {
      assert(N == size() && "compatible sizes");
      array_base<E, N, S> out;
      std::copy(begin(), end(), out.begin());
      return out;
    }
  };
  template <class T>
  std::ostream &operator<<(std::ostream &os,

                           types::dynamic_tuple<T> const &v)
  {
    os << '(';
    size_t n = v.size();
    if (n) {
      os << v.fast(0);
      for (size_t i = 1; i < n; ++i)
        os << ", " << v.fast(i);
    }
    return os << ')';
  }
} // namespace types

PYTHONIC_NS_END

namespace std
{

  template <size_t I, class T>
  typename pythonic::types::dynamic_tuple<T>::const_reference
  get(pythonic::types::dynamic_tuple<T> const &t)
  {
    return t[I];
  }

  template <size_t I, class T>
  typename pythonic::types::dynamic_tuple<T>::reference get(pythonic::types::dynamic_tuple<T> &t)
  {
    return t[I];
  }

  template <size_t I, class T>
  typename pythonic::types::dynamic_tuple<T>::reference get(pythonic::types::dynamic_tuple<T> &&t)
  {
    return t[I];
  }

  template <size_t I, class T>
  struct tuple_element<I, pythonic::types::dynamic_tuple<T>> {
    using type = typename pythonic::types::dynamic_tuple<T>::value_type;
  };
} // namespace std

/* specialize std::hash */
namespace std
{
  template <class T>
  struct hash<pythonic::types::dynamic_tuple<T>> {
    size_t operator()(pythonic::types::dynamic_tuple<T> const &l) const;
  };
} // namespace std

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/include/utils/fwd.hpp"
#include "pythonic/include/utils/seq.hpp"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <typename T>
struct to_python<types::dynamic_tuple<T>> {
  static PyObject *convert(types::dynamic_tuple<T> const &t);
};

PYTHONIC_NS_END
#endif

#endif
