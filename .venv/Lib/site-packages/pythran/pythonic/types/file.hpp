#ifndef PYTHONIC_TYPES_FILE_HPP
#define PYTHONIC_TYPES_FILE_HPP

#include "pythonic/include/types/file.hpp"

#include "pythonic/builtins/IOError.hpp"
#include "pythonic/builtins/RuntimeError.hpp"
#include "pythonic/builtins/StopIteration.hpp"
#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/types/NoneType.hpp"
#include "pythonic/types/assignable.hpp"
#include "pythonic/types/attr.hpp"
#include "pythonic/types/list.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/shared_ref.hpp"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <string>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

PYTHONIC_NS_BEGIN

namespace types
{

  /// _file implementation

  inline _file::_file() : f(nullptr)
  {
  }

  // TODO : no check on file existance?
  inline _file::_file(types::str const &filename, types::str const &strmode)
      : f(fopen(filename.c_str(), strmode.c_str()))
  {
  }

  inline FILE *_file::operator*() const
  {
    return f;
  }

  inline _file::~_file()
  {
    if (f)
      fclose(f);
  }

  /// file implementation

  // Constructors
  inline file::file() : file_iterator(), data(utils::no_memory())
  {
  }

  inline file::file(types::str const &filename, types::str const &strmode)
      : file_iterator(), data(utils::no_memory()), mode(strmode), name(filename), newlines('\n')
  {
    open(filename, strmode);
    if (mode.find_first_of("r+") != -1)
      *(file_iterator *)this = file_iterator(*this);
  }

  // Iterators
  inline file::iterator file::begin()
  {
    return *this;
  }

  inline file::iterator file::end()
  {
    return {};
  }

  // Modifiers
  inline void file::open(types::str const &filename, types::str const &strmode)
  {
    const char *smode = strmode.c_str();
    // Python enforces that the mode, after stripping 'U', begins with 'r',
    // 'w' || 'a'.
    if (*smode == 'U') {
      ++smode;
    } // Not implemented yet

    data = utils::shared_ref<container_type>(filename, smode);
    if (!**data)
      throw types::IOError("Couldn't open file " + filename);
    is_open = true;
  }

  inline void file::close()
  {
    fclose(**data);
    data->f = nullptr;
    is_open = false;
  }

  inline bool file::closed() const
  {
    return !is_open;
  }

  inline types::str const &file::getmode() const
  {
    return mode;
  }

  inline types::str const &file::getname() const
  {
    return name;
  }

  inline types::str const &file::getnewlines() const
  {
    // Python seems to always return none... Doing the same here
    return newlines;
  }

  inline bool file::eof()
  {
    return ::feof(**data);
  }

  inline void file::flush()
  {
    if (!is_open)
      throw ValueError("I/O operation on closed file");
    fflush(**data);
  }

  inline long file::fileno() const
  {
    if (!is_open)
      throw ValueError("I/O operation on closed file");
    return ::fileno(**data);
  }

  inline bool file::isatty() const
  {
    if (!is_open)
      throw ValueError("I/O operation on closed file");
    return ::isatty(this->fileno());
  }

  inline types::str file::read(long size)
  {
    if (!is_open)
      throw ValueError("I/O operation on closed file");
    if (mode.find_first_of("r+") == -1)
      throw IOError("File not open for reading");
    if (size == 0 || (feof(**data) && mode.find_first_of("ra") == -1))
      return types::str();
    long curr_pos = tell();
    seek(0, SEEK_END);
    size = size < 0 ? tell() - curr_pos : size;
    seek(curr_pos);
    char *content = utils::allocate<char>(size);

    // This part needs a new implementation of types::str(char*) to avoid
    // unnecessary copy.
    types::str res(content, fread(content, sizeof(char), size, **data));
    utils::deallocate(content);
    return res;
  }

  template <class T>
  inline void file::read_as(long n, T *buffer)
  {
    if (fread(buffer, sizeof(T), n, **data) < size_t(n)) {
      throw EOFError("read() didn't return enough bytes");
    }
  }

  inline types::str file::readline(long size)
  {
    if (!is_open)
      throw ValueError("I/O operation on closed file");
    if (mode.find_first_of("r+") == -1)
      throw IOError("File not open for reading");
    constexpr static long BUFFER_SIZE = 1024;
    types::str res;
    char read_str[BUFFER_SIZE];

    for (long i = 0; i < size; i += BUFFER_SIZE) {
      // +1 because we read the last chunk so we don't want to count \0
      if (fgets(read_str, std::min(BUFFER_SIZE - 1, size - i) + 1, **data))
        res += read_str;
      if (feof(**data) || res[res.size() - 1] == "\n")
        break;
    }
    return res;
  }

  inline types::list<types::str> file::readlines(long sizehint)
  {
    // Official python doc specifies that sizehint is used as a max of chars
    // But it has not been implemented in the standard python interpreter...
    types::str str;
    types::list<types::str> lst(0);
    while ((str = readline()))
      lst.push_back(str);
    return lst;
  }

  inline void file::seek(long offset, long whence)
  {
    if (!is_open)
      throw ValueError("I/O operation on closed file");
    if (whence != SEEK_SET && whence != SEEK_CUR && whence != SEEK_END)
      throw IOError("file.seek() :  Invalid argument.");
    fseek(**data, offset, whence);
  }

  inline long file::tell() const
  {
    if (!is_open)
      throw ValueError("I/O operation on closed file");
    return ftell(**data);
  }

  inline void file::truncate(long size)
  {
    if (!is_open)
      throw ValueError("I/O operation on closed file");
    if (mode.find_first_of("wa+") == -1)
      throw IOError("file.write() :  File not open for writing.");
    if (size < 0)
      size = this->tell();
#ifdef _WIN32
    long error = _chsize_s(fileno(), size);
#else
    long error = ftruncate(fileno(), size);
#endif
    if (error == -1)
      throw RuntimeError(strerror(errno));
  }

  inline long file::write(types::str const &str)
  {
    if (!is_open)
      throw ValueError("I/O operation on closed file");
    if (mode.find_first_of("wa+") == -1)
      throw IOError("file.write() :  File not open for writing.");
    return fwrite(str.c_str(), sizeof(char), str.size(), **data);
  }

  template <class T>
  void file::writelines(T const &seq)
  {
    auto end = seq.end();
    for (auto it = seq.begin(); it != end; ++it)
      write(*it);
  }

  /// file_iterator implementation
  // TODO : What if the file disappears before the end?
  // Like in :
  // for line in open("myfile"):
  //     print line
  inline file_iterator::file_iterator(file &ref) : f(&ref), set(false), curr(), position(ref.tell())
  {
  }

  inline file_iterator::file_iterator()
      : f(nullptr), set(false), curr(), position(std::numeric_limits<long>::max()) {};

  inline bool file_iterator::operator==(file_iterator const &f2) const
  {
    return position == f2.position;
  }

  inline bool file_iterator::operator!=(file_iterator const &f2) const
  {
    return position != f2.position;
  }

  inline bool file_iterator::operator<(file_iterator const &f2) const
  {
    // Not really elegant...
    // Equivalent to 'return *this != f2;'
    return position < f2.position;
  }

  inline file_iterator &file_iterator::operator++()
  {
    if (f->eof())
      return *this;
    operator*();
    set = false;
    operator*();
    position = f->eof() ? std::numeric_limits<long>::max() : f->tell();
    return *this;
  }

  inline types::str file_iterator::operator*() const
  {
    if (!set) {
      curr = f->readline();
      set = true;
    }
    return curr.chars(); // to make a copy
  }
} // namespace types
PYTHONIC_NS_END

/* pythran attribute system { */
PYTHONIC_NS_BEGIN

namespace builtins
{
  inline bool getattr(types::attr::CLOSED, types::file const &f)
  {
    return f.closed();
  }

  inline types::str const &getattr(types::attr::MODE, types::file const &f)
  {
    return f.getmode();
  }

  inline types::str const &getattr(types::attr::NAME, types::file const &f)
  {
    return f.getname();
  }

  // Python seems to always return none... Doing the same.
  inline types::none_type getattr(types::attr::NEWLINES, types::file const &f)
  {
    return builtins::None;
  }
} // namespace builtins
PYTHONIC_NS_END

/* } */

#endif
