#ifndef PYTHONIC_INCLUDE_TYPES_STR_HPP
#define PYTHONIC_INCLUDE_TYPES_STR_HPP

#include "pythonic/include/types/slice.hpp"
#include "pythonic/include/types/tuple.hpp"

#include "pythonic/include/types/assignable.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/int_.hpp"
#include "pythonic/include/utils/shared_ref.hpp"

#include <cassert>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>

PYTHONIC_NS_BEGIN

namespace types
{

  class str;
  struct const_sliced_str_iterator;

  struct chr {
    char c;
    chr() = default;
    chr(char c) : c(c)
    {
    }
    bool operator==(chr other) const
    {
      return c == other.c;
    }
    long size() const
    {
      return 1;
    }
    operator long() const
    {
      return c - '0';
    }
    operator str() const;
    char const *begin() const
    {
      return &c;
    }
    char const *end() const
    {
      return 1 + &c;
    }
    std::array<char, 1> chars() const
    {
      return {{c}};
    }
  };

  template <class S = slice>
  class sliced_str
  {

    using container_type = std::string;
    utils::shared_ref<container_type> data;

    typename S::normalized_type slicing;

  public:
    //  types
    using reference = container_type::reference;
    using const_reference = container_type::const_reference;
    using iterator = const_sliced_str_iterator;
    using const_iterator = const_sliced_str_iterator;
    using size_type = container_type::size_type;
    using difference_type = container_type::difference_type;
    using value_type = container_type::value_type;
    using allocator_type = container_type::allocator_type;
    using pointer = container_type::pointer;
    using const_pointer = container_type::const_pointer;

    // constructor
    sliced_str();
    sliced_str(sliced_str const &s);

    sliced_str(sliced_str const &s, typename S::normalized_type const &sl);
    sliced_str(types::str const &other, typename S::normalized_type const &s);

    // const getter
    container_type const &get_data() const;
    typename S::normalized_type const &get_slice() const;

    // assignment
    sliced_str &operator=(str const &);
    sliced_str &operator=(sliced_str const &);
    str operator+(sliced_str const &);

    // iterators
    const_iterator begin() const;
    const_iterator end() const;

    // size
    long size() const;

    // accessor
    chr operator[](long i) const;
    chr fast(long i) const;
    template <class Sp>
    std::enable_if_t<is_slice<Sp>::value, sliced_str<Sp>> operator[](Sp const &s) const;

    // conversion
    operator long() const;
    explicit operator bool() const;
    bool operator!() const;

    size_t find(str const &s, size_t pos = 0) const;
    bool contains(str const &v) const;
    bool operator==(str const &v) const;

    // io
    template <class SS>
    friend std::ostream &operator<<(std::ostream &os, types::sliced_str<SS> const &v);
  };

  struct string_iterator;

  class str
  {

    template <class S>
    friend class sliced_str;

    using container_type = std::string;
    utils::shared_ref<container_type> data;

  public:
    static const size_t npos = -1 /*std::string::npos*/;
    static constexpr bool is_vectorizable = false;

    using value_type = str; // in Python, a string contains... strings
    using iterator = string_iterator;
    using reverse_iterator = std::reverse_iterator<string_iterator>;
    using const_reverse_iterator = std::reverse_iterator<string_iterator>;

    str();
    str(std::string const &s);
    str(std::string &&s);
    explicit str(char c);
    str(const char *s);
    template <size_t N>
    str(const char (&s)[N]);
    str(const char *s, size_t n);
    template <class S>
    str(sliced_str<S> const &other);
    template <class T>
    str(T const &begin, T const &end);
    template <class T>
    explicit str(T const &);

    explicit operator char() const;
    explicit operator long int() const;
    explicit operator float() const;
    explicit operator double() const;

    template <class S>
    str &operator=(sliced_str<S> const &other);

    types::str &operator+=(types::str const &s);
    types::str &operator+=(types::chr const &s);

    long size() const;
    iterator begin() const;
    reverse_iterator rbegin() const;
    iterator end() const;
    reverse_iterator rend() const;
    auto c_str() const -> decltype(data->c_str());
    container_type &chars()
    {
      return *data;
    }
    container_type const &chars() const
    {
      return *data;
    }
    auto resize(long n) -> decltype(data->resize(n));
    long find(str const &s, size_t pos = 0) const;
    bool contains(str const &v) const;
    long find_first_of(str const &s, size_t pos = 0) const;
    long find_first_of(const char *s, size_t pos = 0) const;
    long find_first_not_of(str const &s, size_t pos = 0) const;
    long find_last_not_of(str const &s, size_t pos = npos) const;
    str substr(size_t pos = 0, size_t len = npos) const;
    bool empty() const;
    int compare(size_t pos, size_t len, str const &str) const;
    void reserve(size_t n);
    str &replace(size_t pos, size_t len, str const &str);

    template <class S>
    str &operator+=(sliced_str<S> const &other);
    bool operator==(str const &other) const;
    bool operator!=(str const &other) const;
    bool operator<=(str const &other) const;
    bool operator<(str const &other) const;
    bool operator>=(str const &other) const;
    bool operator>(str const &other) const;
    template <class S>
    bool operator==(sliced_str<S> const &other) const;
    bool operator==(chr other) const;

    template <class S>
    std::enable_if_t<is_slice<S>::value, sliced_str<S>> operator()(S const &s) const;

    chr operator[](long i) const;
    chr fast(long i) const;

    template <class S>
    std::enable_if_t<is_slice<S>::value, sliced_str<S>> operator[](S const &s) const;

    explicit operator bool() const;
    long count(types::str const &sub) const;

    intptr_t id() const
    {
      return reinterpret_cast<intptr_t>(&(*data));
    }
  };

  struct string_iterator
      : std::iterator<std::random_access_iterator_tag, str, std::ptrdiff_t, str *, str> {
    std::string::const_iterator curr;
    string_iterator() = default;
    string_iterator(std::string::const_iterator iter) : curr(iter)
    {
    }
    chr operator*() const
    {
      return chr(*curr);
    }
    string_iterator &operator++()
    {
      ++curr;
      return *this;
    }
    string_iterator &operator+=(std::size_t n)
    {
      curr += n;
      return *this;
    }
    string_iterator operator+(std::size_t n)
    {
      return {curr + n};
    }
    string_iterator &operator--()
    {
      --curr;
      return *this;
    }
    string_iterator &operator-=(std::size_t n)
    {
      curr -= n;
      return *this;
    }
    string_iterator operator-(std::size_t n)
    {
      return {curr - n};
    }
    bool operator==(string_iterator const &other) const
    {
      return curr == other.curr;
    }
    bool operator!=(string_iterator const &other) const
    {
      return curr != other.curr;
    }
    bool operator<(string_iterator const &other) const
    {
      return curr < other.curr;
    }
    bool operator<=(string_iterator const &other) const
    {
      return curr <= other.curr;
    }
    std::ptrdiff_t operator-(string_iterator const &other) const
    {
      return curr - other.curr;
    }
  };
  struct const_sliced_str_iterator
      : std::iterator<std::random_access_iterator_tag, str, std::ptrdiff_t, str *, str> {
    const char *data;
    long step;
    const_sliced_str_iterator(char const *data, long step);
    const_sliced_str_iterator operator++();
    bool operator<(const_sliced_str_iterator const &other) const;
    bool operator==(const_sliced_str_iterator const &other) const;
    bool operator!=(const_sliced_str_iterator const &other) const;
    chr operator*() const;
    const_sliced_str_iterator operator-(long n) const;
    long operator-(const_sliced_str_iterator const &other) const;
  };

  str operator+(str const &self, str const &other);
  str operator+(chr const &self, chr const &other);
  str operator+(chr const &self, str const &other);
  str operator+(chr const &self, str const &other);

  template <size_t N>
  str operator+(str const &self, char const (&other)[N]);
  template <size_t N>
  str operator+(chr const &self, char const (&other)[N]);

  template <size_t N>
  str operator+(char const (&self)[N], str const &other);
  template <size_t N>
  str operator+(char const (&self)[N], chr const &other);

  template <size_t N>
  bool operator==(char const (&self)[N], str const &other);

  bool operator==(chr self, str const &other);

  std::ostream &operator<<(std::ostream &os, chr const &s);
  std::ostream &operator<<(std::ostream &os, str const &s);

  str operator*(str const &s, long n);
  str operator*(long t, str const &s);
  str operator*(chr const &s, long n);
  str operator*(long t, chr const &s);
} // namespace types

namespace operator_
{

  template <size_t N, class Arg>
  auto mod(const char (&fmt)[N], Arg &&arg)
      -> decltype(pythonic::types::str(fmt) % std::forward<Arg>(arg));

  pythonic::types::str add(char const *self, char const *other);

  pythonic::types::str mul(char const *self, long other);
  pythonic::types::str mul(long self, char const *other);
} // namespace operator_

template <>
struct assignable<types::chr> {
  using type = types::str;
};

template <>
struct assignable<char *> {
  using type = types::str;
};
template <>
struct assignable<char const *> {
  using type = types::str;
};
template <size_t N>
struct assignable<char[N]> {
  using type = types::str;
};
template <size_t N>
struct assignable<char const[N]> {
  using type = types::str;
};
PYTHONIC_NS_END

namespace std
{
  template <>
  struct hash<pythonic::types::str> {
    size_t operator()(const pythonic::types::str &x) const;
  };

  template <>
  struct hash<pythonic::types::chr> {
    size_t operator()(const pythonic::types::chr &x) const;
  };

  /* std::get overload */
  template <size_t I>
  pythonic::types::str get(pythonic::types::str const &t);

  template <size_t I>
  struct tuple_element<I, pythonic::types::str> {
    using type = pythonic::types::str;
  };

  template <size_t I, class S>
  struct tuple_element<I, pythonic::types::sliced_str<S>> {
    using type = pythonic::types::str;
  };
} // namespace std

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <>
struct __combined<char const *, pythonic::types::str> {
  using type = pythonic::types::str;
};

template <>
struct __combined<pythonic::types::str, char const *> {
  using type = pythonic::types::str;
};

template <size_t N, size_t M>
struct __combined<char[N], char[M]> {
  using type = pythonic::types::str;
};

/* } */
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<types::str> {
  static PyObject *convert(types::str const &v);
};

template <>
struct to_python<types::chr> {
  static PyObject *convert(types::chr const &v);
};

template <class S>
struct to_python<types::sliced_str<S>> {
  static PyObject *convert(types::sliced_str<S> const &v);
};

template <>
struct from_python<types::str> {
  static bool is_convertible(PyObject *obj);
  static types::str convert(PyObject *obj);
};
PYTHONIC_NS_END

#endif

#endif
