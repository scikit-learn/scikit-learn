#ifndef PYTHONIC_INCLUDE_TYPES_FILE_HPP
#define PYTHONIC_INCLUDE_TYPES_FILE_HPP

#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/types/assignable.hpp"
#include "pythonic/include/types/attr.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/shared_ref.hpp"

#include <fstream>
#include <iterator>
#include <string>

PYTHONIC_NS_BEGIN

namespace types
{
  class file;

  struct file_iterator : std::iterator<std::forward_iterator_tag, types::str, ptrdiff_t,
                                       types::str *, types::str /* no ref */> {
  private:
    file *f;
    mutable bool set;
    mutable types::str curr;
    long position;

  public:
    using value_type = types::str;

    file_iterator(file &ref);
    file_iterator();
    bool operator==(file_iterator const &f2) const;
    bool operator!=(file_iterator const &f2) const;
    bool operator<(file_iterator const &f2) const;
    file_iterator &operator++();
    types::str operator*() const;
  };

  struct _file {
    FILE *f;
    _file();
    _file(types::str const &filename, types::str const &strmode = "r");
    FILE *operator*() const;
    ~_file();
  };

  class file : public file_iterator
  {

  private:
    using container_type = _file;
    utils::shared_ref<container_type> data;
    bool is_open;
    types::str mode, name, newlines;

  public:
    // Types
    using iterator = file_iterator;
    using value_type = types::str;

    // Constructors
    file();
    file(types::str const &filename, types::str const &strmode = "r");

    // Iterators
    iterator begin();
    iterator end();

    // Modifiers
    void open(types::str const &filename, types::str const &strmode);

    void close();

    bool closed() const;

    types::str const &getmode() const;

    types::str const &getname() const;

    types::str const &getnewlines() const;

    bool eof();

    void flush();

    long fileno() const;

    bool isatty() const;

    types::str read(long size = -1);

    template <class T>
    void read_as(long n, T *buffer);

    types::str readline(long size = std::numeric_limits<long>::max());

    types::list<types::str> readlines(long sizehint = -1);

    void seek(long offset, long whence = SEEK_SET);

    long tell() const;

    void truncate(long size = -1);

    long write(types::str const &str);

    template <class T>
    void writelines(T const &seq);
  };
} // namespace types
PYTHONIC_NS_END

/* pythran attribute system { */
PYTHONIC_NS_BEGIN

namespace builtins
{
  bool getattr(types::attr::CLOSED, types::file const &f);

  types::str const &getattr(types::attr::MODE, types::file const &f);

  types::str const &getattr(types::attr::NAME, types::file const &f);

  // Python seems to always return none... Doing the same.
  types::none_type getattr(types::attr::NEWLINES, types::file const &f);
} // namespace builtins
PYTHONIC_NS_END

/* } */

#endif
