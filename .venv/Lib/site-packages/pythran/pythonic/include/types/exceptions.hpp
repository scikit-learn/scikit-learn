#ifndef PYTHONIC_INCLUDE_TYPES_EXCEPTIONS_HPP
#define PYTHONIC_INCLUDE_TYPES_EXCEPTIONS_HPP

#include "pythonic/include/builtins/str.hpp"
#include "pythonic/include/types/attr.hpp"
#include "pythonic/include/types/dynamic_tuple.hpp"
#include "pythonic/include/types/str.hpp"

#include <stdexcept>

PYTHONIC_NS_BEGIN

namespace types
{

  class BaseException : public std::exception
  {
  public:
    BaseException(const BaseException &e) = default;
    template <typename... Types>
    BaseException(Types const &...types);
    virtual ~BaseException() noexcept = default;
    dynamic_tuple<str> args;
  };

// Use this to create a python exception class
#define CLASS_EXCEPTION_DECL(name, parent)                                                         \
  class name : public parent                                                                       \
  {                                                                                                \
  public:                                                                                          \
    name() = default;                                                                              \
    name(const name &e) = default;                                                                 \
    template <class... Types>                                                                      \
    name(Types const &...types) : parent(types...)                                                 \
    {                                                                                              \
    }                                                                                              \
    virtual ~name() noexcept = default;                                                            \
  };

  CLASS_EXCEPTION_DECL(SystemExit, BaseException);
  CLASS_EXCEPTION_DECL(KeyboardInterrupt, BaseException);
  CLASS_EXCEPTION_DECL(GeneratorExit, BaseException);
  CLASS_EXCEPTION_DECL(Exception, BaseException);
  CLASS_EXCEPTION_DECL(StopIteration, Exception);
  CLASS_EXCEPTION_DECL(StandardError, Exception);
  CLASS_EXCEPTION_DECL(Warning, Exception);
  CLASS_EXCEPTION_DECL(BytesWarning, Warning);
  CLASS_EXCEPTION_DECL(UnicodeWarning, Warning);
  CLASS_EXCEPTION_DECL(ImportWarning, Warning);
  CLASS_EXCEPTION_DECL(FutureWarning, Warning);
  CLASS_EXCEPTION_DECL(UserWarning, Warning);
  CLASS_EXCEPTION_DECL(SyntaxWarning, Warning);
  CLASS_EXCEPTION_DECL(RuntimeWarning, Warning);
  CLASS_EXCEPTION_DECL(PendingDeprecationWarning, Warning);
  CLASS_EXCEPTION_DECL(DeprecationWarning, Warning);
  CLASS_EXCEPTION_DECL(BufferError, StandardError);
  CLASS_EXCEPTION_DECL(FileNotFoundError, StandardError);
  CLASS_EXCEPTION_DECL(ArithmeticError, StandardError);
  CLASS_EXCEPTION_DECL(AssertionError, StandardError);
  CLASS_EXCEPTION_DECL(AttributeError, StandardError);
  CLASS_EXCEPTION_DECL(EnvironmentError, StandardError);
  CLASS_EXCEPTION_DECL(EOFError, StandardError);
  CLASS_EXCEPTION_DECL(ImportError, StandardError);
  CLASS_EXCEPTION_DECL(LookupError, StandardError);
  CLASS_EXCEPTION_DECL(MemoryError, StandardError);
  CLASS_EXCEPTION_DECL(NameError, StandardError);
  CLASS_EXCEPTION_DECL(ReferenceError, StandardError);
  CLASS_EXCEPTION_DECL(RuntimeError, StandardError);
  CLASS_EXCEPTION_DECL(SyntaxError, StandardError);
  CLASS_EXCEPTION_DECL(SystemError, StandardError);
  CLASS_EXCEPTION_DECL(TypeError, StandardError);
  CLASS_EXCEPTION_DECL(ValueError, StandardError);
  CLASS_EXCEPTION_DECL(FloatingPointError, ArithmeticError);
  CLASS_EXCEPTION_DECL(OverflowError, ArithmeticError);
  CLASS_EXCEPTION_DECL(ZeroDivisionError, ArithmeticError);
  CLASS_EXCEPTION_DECL(IOError, EnvironmentError);
  CLASS_EXCEPTION_DECL(OSError, EnvironmentError);
  CLASS_EXCEPTION_DECL(WindowsError, OSError);
  CLASS_EXCEPTION_DECL(VMSError, OSError);
  CLASS_EXCEPTION_DECL(IndexError, LookupError);
  CLASS_EXCEPTION_DECL(KeyError, LookupError);
  CLASS_EXCEPTION_DECL(UnboundLocalError, NameError);
  CLASS_EXCEPTION_DECL(NotImplementedError, RuntimeError);
  CLASS_EXCEPTION_DECL(IndentationError, SyntaxError);
  CLASS_EXCEPTION_DECL(TabError, IndentationError);
  CLASS_EXCEPTION_DECL(UnicodeError, ValueError);
} // namespace types
PYTHONIC_NS_END

#include "pythonic/include/utils/functor.hpp"
#define PYTHONIC_EXCEPTION_DECL(name)                                                              \
  template <typename... Types>                                                                     \
  types::name name(Types const &...args);                                                          \
                                                                                                   \
  DEFINE_FUNCTOR(pythonic::builtins, name);

/* pythran attribute system { */
#define DECLARE_EXCEPTION_GETATTR(name)                                                            \
  PYTHONIC_NS_BEGIN                                                                                \
  namespace builtins                                                                               \
  {                                                                                                \
    types::none<types::dynamic_tuple<types::str>> getattr(types::attr::ARGS,                       \
                                                          types::name const &f);                   \
  }                                                                                                \
  PYTHONIC_NS_END

#define DECLARE_EXCEPTION_GETATTR_FULL(name)                                                       \
  PYTHONIC_NS_BEGIN                                                                                \
  namespace builtins                                                                               \
  {                                                                                                \
    types::none<types::dynamic_tuple<types::str>> getattr(types::attr::ARGS,                       \
                                                          types::name const &e);                   \
    types::none<types::str> getattr(types::attr::ERRNO, types::name const &e);                     \
    types::none<types::str> getattr(types::attr::STRERROR, types::name const &e);                  \
    types::none<types::str> getattr(types::attr::FILENAME, types::name const &e);                  \
  }                                                                                                \
  PYTHONIC_NS_END

DECLARE_EXCEPTION_GETATTR(BaseException);
DECLARE_EXCEPTION_GETATTR(SystemExit);
DECLARE_EXCEPTION_GETATTR(KeyboardInterrupt);
DECLARE_EXCEPTION_GETATTR(GeneratorExit);
DECLARE_EXCEPTION_GETATTR(Exception);
DECLARE_EXCEPTION_GETATTR(StopIteration);
DECLARE_EXCEPTION_GETATTR(StandardError);
DECLARE_EXCEPTION_GETATTR(Warning);
DECLARE_EXCEPTION_GETATTR(BytesWarning);
DECLARE_EXCEPTION_GETATTR(UnicodeWarning);
DECLARE_EXCEPTION_GETATTR(ImportWarning);
DECLARE_EXCEPTION_GETATTR(FutureWarning);
DECLARE_EXCEPTION_GETATTR(UserWarning);
DECLARE_EXCEPTION_GETATTR(SyntaxWarning);
DECLARE_EXCEPTION_GETATTR(RuntimeWarning);
DECLARE_EXCEPTION_GETATTR(PendingDeprecationWarning);
DECLARE_EXCEPTION_GETATTR(DeprecationWarning);
DECLARE_EXCEPTION_GETATTR(BufferError);
DECLARE_EXCEPTION_GETATTR(FileNotFoundError);
DECLARE_EXCEPTION_GETATTR(ArithmeticError);
DECLARE_EXCEPTION_GETATTR(AssertionError);
DECLARE_EXCEPTION_GETATTR(AttributeError);
DECLARE_EXCEPTION_GETATTR(EOFError);
DECLARE_EXCEPTION_GETATTR(ImportError);
DECLARE_EXCEPTION_GETATTR(LookupError);
DECLARE_EXCEPTION_GETATTR(MemoryError);
DECLARE_EXCEPTION_GETATTR(NameError);
DECLARE_EXCEPTION_GETATTR(ReferenceError);
DECLARE_EXCEPTION_GETATTR(RuntimeError);
DECLARE_EXCEPTION_GETATTR(SyntaxError);
DECLARE_EXCEPTION_GETATTR(SystemError);
DECLARE_EXCEPTION_GETATTR(TypeError);
DECLARE_EXCEPTION_GETATTR(ValueError);
DECLARE_EXCEPTION_GETATTR(FloatingPointError);
DECLARE_EXCEPTION_GETATTR(OverflowError);
DECLARE_EXCEPTION_GETATTR(ZeroDivisionError);
DECLARE_EXCEPTION_GETATTR(IndexError);
DECLARE_EXCEPTION_GETATTR(KeyError);
DECLARE_EXCEPTION_GETATTR(UnboundLocalError);
DECLARE_EXCEPTION_GETATTR(NotImplementedError);
DECLARE_EXCEPTION_GETATTR(IndentationError);
DECLARE_EXCEPTION_GETATTR(TabError);
DECLARE_EXCEPTION_GETATTR(UnicodeError);
DECLARE_EXCEPTION_GETATTR_FULL(IOError);
DECLARE_EXCEPTION_GETATTR_FULL(EnvironmentError);
DECLARE_EXCEPTION_GETATTR_FULL(OSError);

PYTHONIC_NS_BEGIN

namespace types
{

  std::ostream &operator<<(std::ostream &o, BaseException const &e);

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
  std::ostream &operator<<(std::ostream &o, EnvironmentError const &e);
} // namespace types
PYTHONIC_NS_END

/* } */

#endif
