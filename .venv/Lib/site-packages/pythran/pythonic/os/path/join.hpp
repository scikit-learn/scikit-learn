#ifndef PYTHONIC_OS_PATH_JOIN_HPP
#define PYTHONIC_OS_PATH_JOIN_HPP

#ifdef WIN32
#define OS_SEP '\\'
#else
#define OS_SEP '/'
#endif

#include "pythonic/include/os/path/join.hpp"

#include "pythonic/types/str.hpp"

PYTHONIC_NS_BEGIN
namespace os
{
  namespace path
  {
    namespace
    {
      template <class T>
      size_t sizeof_string(T const &s)
      {
        return s.size();
      }

      template <class T, class... Types>
      size_t sizeof_string(T const &s, Types &&...tail)
      {
        return s.size() + sizeof_string(std::forward<Types>(tail)...);
      }

      inline void _join(types::str &buffer)
      {
      }

      template <class T, class... Types>
      void _join(types::str &buffer, T &&head, Types &&...tail)
      {
        if (((types::str)head)[0] == "/")
          buffer = std::forward<T>(head);
        else if (!buffer || *buffer.chars().rbegin() == OS_SEP || *buffer.rbegin() == "/")
          buffer += std::forward<T>(head);
        else {
          buffer.chars() += OS_SEP;
          buffer += std::forward<T>(head);
        }
        _join(buffer, std::forward<Types>(tail)...);
      }
    } // namespace

    template <class T>
    T join(T &&head)
    {
      return head;
    }

    template <class T, class... Types>
    types::str join(T &&head, Types &&...tail)
    {
      types::str p = head;
      p.reserve(sizeof_string(tail...));
      _join(p, std::forward<Types>(tail)...);
      return p;
    }
  } // namespace path
} // namespace os
PYTHONIC_NS_END

#endif
