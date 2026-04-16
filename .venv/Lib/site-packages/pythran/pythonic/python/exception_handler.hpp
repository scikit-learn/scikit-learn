#ifndef PYTHONIC_PYTHON_EXCEPTION_HANDLER_HPP
#define PYTHONIC_PYTHON_EXCEPTION_HANDLER_HPP

#ifdef ENABLE_PYTHON_MODULE

#include "Python.h"
#include <utility>

PYTHONIC_NS_BEGIN

// This function have to be include after every others exceptions to have
// correct exception macro defined.
template <class F>
PyObject *handle_python_exception(F &&f)
{
  try {
    return f();
  }
#ifdef PYTHONIC_BUILTIN_SYNTAXWARNING_HPP
  catch (pythonic::types::SyntaxWarning &e) {
    PyErr_SetString(PyExc_SyntaxWarning, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_RUNTIMEWARNING_HPP
  catch (pythonic::types::RuntimeWarning &e) {
    PyErr_SetString(PyExc_RuntimeWarning, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_DEPRECATIONWARNING_HPP
  catch (pythonic::types::DeprecationWarning &e) {
    PyErr_SetString(PyExc_DeprecationWarning, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_IMPORTWARNING_HPP
  catch (pythonic::types::ImportWarning &e) {
    PyErr_SetString(PyExc_ImportWarning, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_UNICODEWARNING_HPP
  catch (pythonic::types::UnicodeWarning &e) {
    PyErr_SetString(PyExc_UnicodeWarning, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_BYTESWARNING_HPP
  catch (pythonic::types::BytesWarning &e) {
    PyErr_SetString(PyExc_BytesWarning, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_USERWARNING_HPP
  catch (pythonic::types::UserWarning &e) {
    PyErr_SetString(PyExc_UserWarning, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_FUTUREWARNING_HPP
  catch (pythonic::types::FutureWarning &e) {
    PyErr_SetString(PyExc_FutureWarning, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_PENDINGDEPRECATIONWARNING_HPP
  catch (pythonic::types::PendingDeprecationWarning &e) {
    PyErr_SetString(PyExc_PendingDeprecationWarning,
                    pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_WARNING_HPP
  catch (pythonic::types::Warning &e) {
    PyErr_SetString(PyExc_Warning, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_UNICODEERROR_HPP
  catch (pythonic::types::UnicodeError &e) {
    PyErr_SetString(PyExc_UnicodeError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_VALUEERROR_HPP
  catch (pythonic::types::ValueError &e) {
    PyErr_SetString(PyExc_ValueError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_TYPEERROR_HPP
  catch (pythonic::types::TypeError &e) {
    PyErr_SetString(PyExc_TypeError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_SYSTEMERROR_HPP
  catch (pythonic::types::SystemError &e) {
    PyErr_SetString(PyExc_SystemError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_TABERROR_HPP
  catch (pythonic::types::TabError &e) {
    PyErr_SetString(PyExc_TabError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_INDENTATIONERROR_HPP
  catch (pythonic::types::IndentationError &e) {
    PyErr_SetString(PyExc_IndentationError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_SYNTAXERROR_HPP
  catch (pythonic::types::SyntaxError &e) {
    PyErr_SetString(PyExc_SyntaxError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_NOTIMPLEMENTEDERROR_HPP
  catch (pythonic::types::NotImplementedError &e) {
    PyErr_SetString(PyExc_NotImplementedError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_RUNTIMEERROR_HPP
  catch (pythonic::types::RuntimeError &e) {
    PyErr_SetString(PyExc_RuntimeError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_REFERENCEERROR_HPP
  catch (pythonic::types::ReferenceError &e) {
    PyErr_SetString(PyExc_ReferenceError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_UNBOUNDLOCALERROR_HPP
  catch (pythonic::types::UnboundLocalError &e) {
    PyErr_SetString(PyExc_UnboundLocalError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_NAMEERROR_HPP
  catch (pythonic::types::NameError &e) {
    PyErr_SetString(PyExc_NameError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_MEMORYERROR_HPP
  catch (pythonic::types::MemoryError &e) {
    PyErr_SetString(PyExc_MemoryError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_KEYERROR_HPP
  catch (pythonic::types::KeyError &e) {
    PyErr_SetString(PyExc_KeyError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_INDEXERROR_HPP
  catch (pythonic::types::IndexError &e) {
    PyErr_SetString(PyExc_IndexError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_LOOKUPERROR_HPP
  catch (pythonic::types::LookupError &e) {
    PyErr_SetString(PyExc_LookupError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_IMPORTERROR_HPP
  catch (pythonic::types::ImportError &e) {
    PyErr_SetString(PyExc_ImportError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_EOFERROR_HPP
  catch (pythonic::types::EOFError &e) {
    PyErr_SetString(PyExc_EOFError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_OSERROR_HPP
  catch (pythonic::types::OSError &e) {
    PyErr_SetString(PyExc_OSError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_IOERROR_HPP
  catch (pythonic::types::IOError &e) {
    PyErr_SetString(PyExc_IOError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_ENVIRONMENTERROR_HPP
  catch (pythonic::types::EnvironmentError &e) {
    PyErr_SetString(PyExc_EnvironmentError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_ATTRIBUTEERROR_HPP
  catch (pythonic::types::AttributeError &e) {
    PyErr_SetString(PyExc_AttributeError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_ASSERTIONERROR_HPP
  catch (pythonic::types::AssertionError &e) {
    PyErr_SetString(PyExc_AssertionError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_ZERODIVISIONERROR_HPP
  catch (pythonic::types::ZeroDivisionError &e) {
    PyErr_SetString(PyExc_ZeroDivisionError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_OVERFLOWERROR_HPP
  catch (pythonic::types::OverflowError &e) {
    PyErr_SetString(PyExc_OverflowError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_FLOATINGPOINTERROR_HPP
  catch (pythonic::types::FloatingPointError &e) {
    PyErr_SetString(PyExc_FloatingPointError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_ARITHMETICERROR_HPP
  catch (pythonic::types::ArithmeticError &e) {
    PyErr_SetString(PyExc_ArithmeticError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_FILENOTFOUNDERROR_HPP
  catch (pythonic::types::FileNotFoundError &e) {
    PyErr_SetString(PyExc_FileNotFoundError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_BUFFERERROR_HPP
  catch (pythonic::types::BufferError &e) {
    PyErr_SetString(PyExc_BufferError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_STANDARDERROR_HPP
  catch (pythonic::types::StandardError &e) {
    PyErr_SetString(PyExc_StandardError, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_STOPITERATION_HPP
  catch (pythonic::types::StopIteration &e) {
    PyErr_SetString(PyExc_StopIteration, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_EXCEPTION_HPP
  catch (pythonic::types::Exception &e) {
    PyErr_SetString(PyExc_Exception, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_GENERATOREXIT_HPP
  catch (pythonic::types::GeneratorExit &e) {
    PyErr_SetString(PyExc_GeneratorExit, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_KEYBOARDINTERRUPT_HPP
  catch (pythonic::types::KeyboardInterrupt &e) {
    PyErr_SetString(PyExc_KeyboardInterrupt, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_SYSTEMEXIT_HPP
  catch (pythonic::types::SystemExit &e) {
    PyErr_SetString(PyExc_SystemExit, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
#ifdef PYTHONIC_BUILTIN_BASEEXCEPTION_HPP
  catch (pythonic::types::BaseException &e) {
    PyErr_SetString(PyExc_BaseException, pythonic::builtins::functor::str{}(e.args).c_str());
  }
#endif
  catch (...) {
    PyErr_SetString(PyExc_RuntimeError, "Something happened on the way to heaven");
  }
  return nullptr;
}
PYTHONIC_NS_END

#endif

#endif
