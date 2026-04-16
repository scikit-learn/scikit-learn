#ifndef PYTHONIC_INCLUDE_TYPES_SLICE_HPP
#define PYTHONIC_INCLUDE_TYPES_SLICE_HPP

#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/types/attr.hpp"

#include <limits>
#include <ostream>

PYTHONIC_NS_BEGIN

namespace types
{
  template <class T>
  class bound
  {
    T value_;

    // use a sentinel to store none<long>, it takes less space
    static constexpr long sentinel = std::numeric_limits<T>::min();

  public:
    bound() = default;
    bound(none_type) : value_(sentinel)
    {
    }
    bound(none<T> v) : value_(v.is_none ? sentinel : (T)v)
    {
    }
    bound(T v) : value_(v)
    {
    }

    operator T() const
    {
      return value_;
    }
    operator none<T>() const
    {
      if (value_ == sentinel)
        return none_type();
      else
        return value_;
    }
    bool is_none() const
    {
      return value_ == sentinel;
    }
  };

  struct slice;
  template <long stride>
  struct cstride_slice;
  template <long stride>
  struct cstride_normalized_slice;
  struct fast_contiguous_slice;

  struct normalized_slice {
    long lower, upper, step;
    normalized_slice() = default;
    normalized_slice(long lower, long upper, long step = 1);

    normalized_slice operator*(normalized_slice const &other) const;
    template <long stride>
    normalized_slice operator*(cstride_normalized_slice<stride> const &other) const;
    normalized_slice operator*(slice const &other) const;
    template <long stride>
    normalized_slice operator*(cstride_slice<stride> const &other) const;
    normalized_slice operator*(fast_contiguous_slice const &other) const;

    long size() const;
    inline long get(long i) const;
  };

  struct slice {
    using normalized_type = normalized_slice;

    bound<long> lower, upper, step;
    slice(none<long> lower, none<long> upper, none<long> step);
    slice();

    slice operator*(slice const &other) const;

    template <long stride>
    slice operator*(cstride_slice<stride> const &other) const;

    slice operator*(fast_contiguous_slice const &other) const;

    /*
       Normalize change a[:-1] to a[:len(a)-1] to have positif index.
       It also check for value bigger than len(a) to fit the size of the
       container
       */
    normalized_slice normalize(long max_size) const;

    /*
     * An assert is raised when we can't compute the size without more
     * informations.
     */
    long size() const;

    long get(long i) const;
  };

  template <long stride>
  struct cstride_normalized_slice {
    long lower, upper;
    static constexpr long step = stride;
    cstride_normalized_slice() = default;
    cstride_normalized_slice(long lower, long upper, long = 0);

    normalized_slice operator*(normalized_slice const &other) const;
    template <long other_stride>
    std::conditional_t<(stride < 256 && other_stride < 256),
                       cstride_normalized_slice<stride * other_stride>, normalized_slice>
    operator*(cstride_normalized_slice<other_stride> const &other) const;

    normalized_slice operator*(slice const &other) const;
    template <long other_stride>
    std::conditional_t<(stride < 256 && other_stride < 256),
                       cstride_normalized_slice<stride * other_stride>, normalized_slice>
    operator*(cstride_slice<other_stride> const &other) const;
    cstride_normalized_slice operator*(fast_contiguous_slice const &other) const;

    long size() const;
    inline long get(long i) const;
  };

  template <long stride>
  constexpr long cstride_normalized_slice<stride>::step;

  template <long stride>
  struct cstride_slice {
    using normalized_type = cstride_normalized_slice<stride>;

    bound<long> lower, upper;
    static constexpr long step = stride;

    cstride_slice(none<long> lower, none<long> upper);
    cstride_slice();

    slice operator*(slice const &other) const;

    template <long other_stride>
    std::conditional_t<(stride < 256 && other_stride < 256), cstride_slice<stride * other_stride>,
                       slice>
    operator*(cstride_slice<other_stride> const &other) const;

    cstride_slice operator*(fast_contiguous_slice const &other) const;

    /*
       Normalize change a[:-1] to a[:len(a)-1] to have positif index.
       It also check for value bigger than len(a) to fit the size of the
       container
       */
    cstride_normalized_slice<stride> normalize(long max_size) const;

    /*
     * An assert is raised when we can't compute the size without more
     * informations.
     */
    long size() const;

    long get(long i) const;
  };

  template <long stride>
  constexpr long cstride_slice<stride>::step;

  // Kept for compatibility with former pythran version -cython use it.
  using contiguous_slice = cstride_slice<1>;

  struct fast_contiguous_slice {
    using normalized_type = cstride_normalized_slice<1>;
    long lower;
    bound<long> upper;
    static constexpr long step = 1;
    fast_contiguous_slice(none<long> lower, none<long> upper);
    fast_contiguous_slice() = default;

    fast_contiguous_slice operator*(fast_contiguous_slice const &other) const;
    template <long stride>
    cstride_slice<stride> operator*(cstride_slice<stride> const &other) const;
    slice operator*(slice const &other) const;

    cstride_normalized_slice<1> normalize(long max_size) const;
    long size() const;
  };

  template <class T>
  struct normalized {
    using type = T;
  };

  template <>
  struct normalized<slice> {
    using type = normalized_slice;
  };

  template <long stride>
  struct normalized<cstride_slice<stride>> {
    using type = cstride_normalized_slice<stride>;
  };

  template <>
  struct normalized<fast_contiguous_slice> {
    using type = cstride_normalized_slice<1>;
  };

  template <class S>
  struct is_slice : std::false_type {
  };
  template <>
  struct is_slice<fast_contiguous_slice> : std::true_type {
  };
  template <>
  struct is_slice<slice> : std::true_type {
  };
  template <long stride>
  struct is_slice<cstride_slice<stride>> : std::true_type {
  };

  template <class S>
  struct is_normalized_slice : std::false_type {
  };
  template <>
  struct is_normalized_slice<normalized_slice> : std::true_type {
  };
  template <long stride>
  struct is_normalized_slice<cstride_normalized_slice<stride>> : std::true_type {
  };

  template <class S>
  using normalize_t = typename normalized<S>::type;

  template <class S>
  std::enable_if_t<!is_slice<S>::value, S> normalize(S s, long n)
  {
    if (s < 0)
      s += n;
    return s;
  }
  inline none_type normalize(none_type s, long n)
  {
    return {};
  }
  template <class S>
  auto normalize(S s, long n) -> decltype(s.normalize(n))
  {
    return s.normalize(n);
  }

  template <class I0, class I1>
  none_type adapt_slice(none_type, I0 const &, I1 const &)
  {
    return {};
  }
  template <class I0, class I1>
  long adapt_slice(long l, I0 const &index0, I1 const &index1)
  {

    if ((long)index0 != (long)index1)
      return 0;
    else
      return l;
  }
  template <class I0, class I1>
  slice adapt_slice(slice const &s, I0 const &index0, I1 const &index1)
  {
    if ((long)index0 != (long)index1)
      return {0, 1, 1};
    else
      return s;
  }
  template <class S, class I0, class I1>
  cstride_slice<1> adapt_slice(S const &s, I0 const &index0, I1 const &index1)
  {
    if ((long)index0 != (long)index1)
      return {0, 1};
    else
      return s;
  }
  template <class I0, class I1>
  fast_contiguous_slice adapt_slice(fast_contiguous_slice const &s, I0 const &index0,
                                    I1 const &index1)
  {
    if ((long)index0 != (long)index1)
      return {0, 1};
    else
      return s;
  }

  template <class S>
  std::enable_if_t<is_slice<S>::value, std::ostream &> operator<<(std::ostream &os, S const &s);
} // namespace types
namespace builtins
{
  template <class T>
  auto getattr(types::attr::START, T const &s) -> decltype(s.lower)
  {
    return s.lower;
  }
  template <class T>
  auto getattr(types::attr::STOP, T const &s) -> decltype(s.upper)
  {
    return s.upper;
  }
  template <class T>
  auto getattr(types::attr::STEP, T const &s) -> decltype(s.step)
  {
    return s.step;
  }
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <class T>
struct to_python<types::bound<T>> {
  static PyObject *convert(types::bound<T> const &n);
};

template <long stride>
struct to_python<types::cstride_slice<stride>> {
  static PyObject *convert(types::cstride_slice<stride> const &n);
};

template <long stride>
struct to_python<types::cstride_normalized_slice<stride>> {
  static PyObject *convert(types::cstride_normalized_slice<stride> const &n);
};

template <>
struct to_python<types::slice> {
  static PyObject *convert(types::slice const &n);
};

template <>
struct to_python<types::normalized_slice> {
  static PyObject *convert(types::normalized_slice const &n);
};

template <>
struct from_python<types::slice> {
  static bool is_convertible(PyObject *obj);
  static types::slice convert(PyObject *obj);
};

PYTHONIC_NS_END

#endif

#endif
