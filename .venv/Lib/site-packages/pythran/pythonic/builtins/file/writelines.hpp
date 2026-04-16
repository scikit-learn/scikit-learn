#ifndef PYTHONIC_BUILTIN_FILE_WRITELINES_HPP
#define PYTHONIC_BUILTIN_FILE_WRITELINES_HPP

#include "pythonic/include/builtins/file/writelines.hpp"

#include "pythonic/types/file.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace file
  {

    template <class F, class T>
    void writelines(F &&f, T const &sequence)
    {
      f.writelines(sequence);
    }
  } // namespace file
} // namespace builtins
PYTHONIC_NS_END
#endif
