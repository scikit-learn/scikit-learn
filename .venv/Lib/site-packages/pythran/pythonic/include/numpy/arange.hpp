#ifndef PYTHONIC_INCLUDE_NUMPY_ARANGE_HPP
#define PYTHONIC_INCLUDE_NUMPY_ARANGE_HPP

#include "pythonic/include/operator_/pos.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {
#ifdef USE_XSIMD
    template <class T>
    struct arange_simd_iterator {
      using vector_type = xsimd::batch<T>;
      vector_type curr_;
      vector_type step_;
      long index_;
      arange_simd_iterator(T start, T step, long n)
          : curr_(), step_(static_cast<T>(vector_type::size * step)),
            index_(static_cast<long>(n / vector_type::size))
      {
        T from[vector_type::size];
        for (size_t i = 0; i < vector_type::size; ++i)
          from[i] = start + i * step;
        curr_ = vector_type::load_unaligned(from);
      }
      vector_type operator*() const
      {
        return curr_;
      }
      arange_simd_iterator &operator++()
      {
        curr_ += step_;
        ++index_;
        return *this;
      }
      arange_simd_iterator &operator+=(long n)
      {
        curr_ += n * step_;
        index_ += n;
        return *this;
      }
      arange_simd_iterator operator+(long n) const
      {
        arange_simd_iterator other{*this};
        return other += n;
      }
      arange_simd_iterator &operator--()
      {
        curr_ -= step_;
        --index_;
        return *this;
      }
      long operator-(arange_simd_iterator const &other) const
      {
        return index_ - other.index_;
      }
      bool operator!=(arange_simd_iterator const &other) const
      {
        return index_ != other.index_;
      }
      bool operator==(arange_simd_iterator const &other) const
      {
        return index_ == other.index_;
      }
      bool operator<(arange_simd_iterator const &other) const
      {
        return index_ < other.index_;
      }
      arange_simd_iterator &operator=(arange_simd_iterator const &other) = default;
    };
#endif
    template <class T>
    struct arange_index {
      T start, step;
      long size;
      using iterator = types::nditerator<arange_index>;
      using const_iterator = types::const_nditerator<arange_index>;
      using dtype = T;
      using value_type = dtype;
      using shape_t = types::pshape<long>;
#ifdef USE_XSIMD
      using simd_iterator = arange_simd_iterator<T>;
      using simd_iterator_nobroadcast = simd_iterator;
      template <class vectorizer>
      simd_iterator vbegin(vectorizer) const
      {
        return {start, step, 0};
      }
      template <class vectorizer>
      simd_iterator vend(vectorizer) const
      {
        return {static_cast<T>(start + size * step), step, size};
      }
#endif
      static constexpr size_t value = 1;
      static constexpr bool is_strided = false;
      static constexpr bool is_vectorizable = types::is_vectorizable<T>::value;

      T fast(long i) const
      {
        return start + i * step;
      }

      dtype load(long i) const
      {
        return fast(i);
      }

      template <size_t I>
      long shape() const
      {
        return size;
      }
      types::ndarray<dtype, shape_t> operator[](types::slice s) const
      {
        auto ns = s.normalize(size);
        arange_index r{start + s.lower * step, step * ns.step, ns.size()};
        return {types::numpy_expr<pythonic::operator_::functor::pos, arange_index>{r}};
      }
      types::ndarray<dtype, shape_t> operator()(types::slice s) const
      {
        return operator[](s);
      }

      template <long stride>
      types::ndarray<dtype, shape_t> operator[](types::cstride_slice<stride> s) const
      {
        auto ns = s.normalize(size);
        arange_index r{start + s.lower * step, step * stride, ns.size()};
        return {types::numpy_expr<pythonic::operator_::functor::pos, arange_index>{r}};
      }

      template <long stride>
      types::ndarray<dtype, shape_t> operator()(types::cstride_slice<stride> s) const
      {
        return operator[](s);
      }

      template <class... S>
      auto operator()(S const &...s) const
          -> std::enable_if_t<(sizeof...(S) > 1),
                              decltype(std::declval<types::ndarray<dtype, shape_t>>()(s...))>
      {
        return types::ndarray<dtype, shape_t>{
            types::numpy_expr<pythonic::operator_::functor::pos, arange_index>{*this}}(s...);
      }
      const_iterator begin() const
      {
        return {*this, 0};
      }
      const_iterator end() const
      {
        return {*this, size};
      }

      iterator begin()
      {
        return {*this, 0};
      }
      iterator end()
      {
        return {*this, size};
      }
    };
  } // namespace details

  template <class T, class U, class S = long,
            class dtype = types::dtype_t<typename __combined<T, U, S>::type>>
  types::numpy_expr<pythonic::operator_::functor::pos, details::arange_index<typename dtype::type>>
  arange(T begin, U end, S step = S(1), dtype d = dtype());

  template <class T>
  types::numpy_expr<pythonic::operator_::functor::pos,
                    details::arange_index<typename types::dtype_t<T>::type>>
  arange(T end);

  DEFINE_FUNCTOR(pythonic::numpy, arange);
} // namespace numpy
PYTHONIC_NS_END

#endif
