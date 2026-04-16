#ifndef PYTHONIC_TYPES_EXCEPTIONS_HPP
#define PYTHONIC_TYPES_EXCEPTIONS_HPP

#include "pythonic/include/types/exceptions.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/builtins/str.hpp"
#include "pythonic/types/attr.hpp"
#include "pythonic/types/dynamic_tuple.hpp"
#include "pythonic/types/str.hpp"

#include <stdexcept>

PYTHONIC_NS_BEGIN

namespace types
{

  template <typename... Types>
  BaseException::BaseException(Types const &...types) : args({builtins::functor::str{}(types)...})
  {
  }

// Use this to create a python exception class
#define CLASS_EXCEPTION_IMPL(name, parent)

  CLASS_EXCEPTION_IMPL(SystemExit, BaseException);
  CLASS_EXCEPTION_IMPL(KeyboardInterrupt, BaseException);
  CLASS_EXCEPTION_IMPL(GeneratorExit, BaseException);
  CLASS_EXCEPTION_IMPL(Exception, BaseException);
  CLASS_EXCEPTION_IMPL(StopIteration, Exception);
  CLASS_EXCEPTION_IMPL(StandardError, Exception);
  CLASS_EXCEPTION_IMPL(Warning, Exception);
  CLASS_EXCEPTION_IMPL(BytesWarning, Warning);
  CLASS_EXCEPTION_IMPL(UnicodeWarning, Warning);
  CLASS_EXCEPTION_IMPL(ImportWarning, Warning);
  CLASS_EXCEPTION_IMPL(FutureWarning, Warning);
  CLASS_EXCEPTION_IMPL(UserWarning, Warning);
  CLASS_EXCEPTION_IMPL(SyntaxWarning, Warning);
  CLASS_EXCEPTION_IMPL(RuntimeWarning, Warning);
  CLASS_EXCEPTION_IMPL(PendingDeprecationWarning, Warning);
  CLASS_EXCEPTION_IMPL(DeprecationWarning, Warning);
  CLASS_EXCEPTION_IMPL(BufferError, StandardError);
  CLASS_EXCEPTION_IMPL(FileNotFoundError, StandardError);
  CLASS_EXCEPTION_IMPL(ArithmeticError, StandardError);
  CLASS_EXCEPTION_IMPL(AssertionError, StandardError);
  CLASS_EXCEPTION_IMPL(AttributeError, StandardError);
  CLASS_EXCEPTION_IMPL(EnvironmentError, StandardError);
  CLASS_EXCEPTION_IMPL(EOFError, StandardError);
  CLASS_EXCEPTION_IMPL(ImportError, StandardError);
  CLASS_EXCEPTION_IMPL(LookupError, StandardError);
  CLASS_EXCEPTION_IMPL(MemoryError, StandardError);
  CLASS_EXCEPTION_IMPL(NameError, StandardError);
  CLASS_EXCEPTION_IMPL(ReferenceError, StandardError);
  CLASS_EXCEPTION_IMPL(RuntimeError, StandardError);
  CLASS_EXCEPTION_IMPL(SyntaxError, StandardError);
  CLASS_EXCEPTION_IMPL(SystemError, StandardError);
  CLASS_EXCEPTION_IMPL(TypeError, StandardError);
  CLASS_EXCEPTION_IMPL(ValueError, StandardError);
  CLASS_EXCEPTION_IMPL(FloatingPointError, ArithmeticError);
  CLASS_EXCEPTION_IMPL(OverflowError, ArithmeticError);
  CLASS_EXCEPTION_IMPL(ZeroDivisionError, ArithmeticError);
  CLASS_EXCEPTION_IMPL(IOError, EnvironmentError);
  CLASS_EXCEPTION_IMPL(OSError, EnvironmentError);
  CLASS_EXCEPTION_IMPL(WindowsError, OSError);
  CLASS_EXCEPTION_IMPL(VMSError, OSError);
  CLASS_EXCEPTION_IMPL(IndexError, LookupError);
  CLASS_EXCEPTION_IMPL(KeyError, LookupError);
  CLASS_EXCEPTION_IMPL(UnboundLocalError, NameError);
  CLASS_EXCEPTION_IMPL(NotImplementedError, RuntimeError);
  CLASS_EXCEPTION_IMPL(IndentationError, SyntaxError);
  CLASS_EXCEPTION_IMPL(TabError, IndentationError);
  CLASS_EXCEPTION_IMPL(UnicodeError, ValueError);
} // namespace types
PYTHONIC_NS_END

#include "pythonic/utils/functor.hpp"
#define PYTHONIC_EXCEPTION_IMPL(name)                                                              \
  template <typename... Types>                                                                     \
  types::name name(Types const &...args)                                                           \
  {                                                                                                \
    return types::name(args...);                                                                   \
  }

/* pythran attribute system { */
#define IMPL_EXCEPTION_GETATTR(name)                                                               \
  PYTHONIC_NS_BEGIN                                                                                \
  namespace builtins                                                                               \
  {                                                                                                \
    inline types::none<types::dynamic_tuple<types::str>> getattr(types::attr::ARGS,                \
                                                                 types::name const &f)             \
    {                                                                                              \
      return f.args;                                                                               \
    }                                                                                              \
  }                                                                                                \
  PYTHONIC_NS_END

#define IMPL_EXCEPTION_GETATTR_FULL(name)                                                          \
  PYTHONIC_NS_BEGIN                                                                                \
  namespace builtins                                                                               \
  {                                                                                                \
    inline types::none<types::dynamic_tuple<types::str>> getattr(types::attr::ARGS,                \
                                                                 types::name const &e)             \
    {                                                                                              \
      if (e.args.size() > 3 || e.args.size() < 2)                                                  \
        return e.args;                                                                             \
      else                                                                                         \
        return types::dynamic_tuple<types::str>(e.args.begin(), e.args.begin() + 2);               \
    }                                                                                              \
    inline types::none<types::str> getattr(types::attr::ERRNO, types::name const &e)               \
    {                                                                                              \
      if (e.args.size() > 3 || e.args.size() < 2)                                                  \
        return builtins::None;                                                                     \
      else                                                                                         \
        return e.args[0];                                                                          \
    }                                                                                              \
    inline types::none<types::str> getattr(types::attr::STRERROR, types::name const &e)            \
    {                                                                                              \
      if (e.args.size() > 3 || e.args.size() < 2)                                                  \
        return builtins::None;                                                                     \
      else                                                                                         \
        return e.args[1];                                                                          \
    }                                                                                              \
    inline types::none<types::str> getattr(types::attr::FILENAME, types::name const &e)            \
    {                                                                                              \
      if (e.args.size() != 3)                                                                      \
        return builtins::None;                                                                     \
      else                                                                                         \
        return e.args[2];                                                                          \
    }                                                                                              \
  }                                                                                                \
  PYTHONIC_NS_END

IMPL_EXCEPTION_GETATTR(BaseException);
IMPL_EXCEPTION_GETATTR(SystemExit);
IMPL_EXCEPTION_GETATTR(KeyboardInterrupt);
IMPL_EXCEPTION_GETATTR(GeneratorExit);
IMPL_EXCEPTION_GETATTR(Exception);
IMPL_EXCEPTION_GETATTR(StopIteration);
IMPL_EXCEPTION_GETATTR(StandardError);
IMPL_EXCEPTION_GETATTR(Warning);
IMPL_EXCEPTION_GETATTR(BytesWarning);
IMPL_EXCEPTION_GETATTR(UnicodeWarning);
IMPL_EXCEPTION_GETATTR(ImportWarning);
IMPL_EXCEPTION_GETATTR(FutureWarning);
IMPL_EXCEPTION_GETATTR(UserWarning);
IMPL_EXCEPTION_GETATTR(SyntaxWarning);
IMPL_EXCEPTION_GETATTR(RuntimeWarning);
IMPL_EXCEPTION_GETATTR(PendingDeprecationWarning);
IMPL_EXCEPTION_GETATTR(DeprecationWarning);
IMPL_EXCEPTION_GETATTR(BufferError);
IMPL_EXCEPTION_GETATTR(FileNotFoundError);
IMPL_EXCEPTION_GETATTR(ArithmeticError);
IMPL_EXCEPTION_GETATTR(AssertionError);
IMPL_EXCEPTION_GETATTR(AttributeError);
IMPL_EXCEPTION_GETATTR(EOFError);
IMPL_EXCEPTION_GETATTR(ImportError);
IMPL_EXCEPTION_GETATTR(LookupError);
IMPL_EXCEPTION_GETATTR(MemoryError);
IMPL_EXCEPTION_GETATTR(NameError);
IMPL_EXCEPTION_GETATTR(ReferenceError);
IMPL_EXCEPTION_GETATTR(RuntimeError);
IMPL_EXCEPTION_GETATTR(SyntaxError);
IMPL_EXCEPTION_GETATTR(SystemError);
IMPL_EXCEPTION_GETATTR(TypeError);
IMPL_EXCEPTION_GETATTR(ValueError);
IMPL_EXCEPTION_GETATTR(FloatingPointError);
IMPL_EXCEPTION_GETATTR(OverflowError);
IMPL_EXCEPTION_GETATTR(ZeroDivisionError);
IMPL_EXCEPTION_GETATTR(IndexError);
IMPL_EXCEPTION_GETATTR(KeyError);
IMPL_EXCEPTION_GETATTR(UnboundLocalError);
IMPL_EXCEPTION_GETATTR(NotImplementedError);
IMPL_EXCEPTION_GETATTR(IndentationError);
IMPL_EXCEPTION_GETATTR(TabError);
IMPL_EXCEPTION_GETATTR(UnicodeError);
IMPL_EXCEPTION_GETATTR_FULL(IOError);
IMPL_EXCEPTION_GETATTR_FULL(EnvironmentError);
IMPL_EXCEPTION_GETATTR_FULL(OSError);

PYTHONIC_NS_BEGIN

namespace types
{

  inline std::ostream &operator<<(std::ostream &o, BaseException const &e)
  {
    return o << e.args;
  }

  /* @brief Convert EnvironmentError to a string.
   *
   * The number of arguments used when creating the EnvironmentError impact
   * the resulting "type" || formatting of the chain. We aim to mimic python
   * behavior of course:
   * - only one arg, then assume it can be converted to string,
   * - two args, then the first one is the errno, the next one a string,
   * - three args, like two args, adding "filename" as third one (after ':')
   * - four || more args, the "tuple" used to construct the exception
   *
   */
  inline std::ostream &operator<<(std::ostream &o, EnvironmentError const &e)
  {
    if (e.args.size() == 1)
      return o << e.args[0];
    if (e.args.size() == 2)
      return o << "[Errno " << e.args[0] << "] " << e.args[1];
    else if (e.args.size() == 3)
      return o << "[Errno " << e.args[0] << "] " << e.args[1] << ": '" << e.args[2] << "'";
    else {
      // Generate "('a', 'b', 'c', 'd') if a,b,c, && d are in e.args
      std::string listsep = "";
      o << "(";
      for (auto &arg : e.args) {
        o << listsep << "'" << arg << "'";
        listsep = ", ";
      }
      o << ")";
      return o;
    }
  }
} // namespace types
PYTHONIC_NS_END

/* } */

#endif
