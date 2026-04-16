#ifndef PYTHONIC_NUMPY_NDARRAY_TOFILE_HPP
#define PYTHONIC_NUMPY_NDARRAY_TOFILE_HPP

#include "pythonic/builtins/FileNotFoundError.hpp"
#include "pythonic/include/numpy/ndarray/tofile.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"
#include <fstream>
#include <limits>

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace ndarray
  {
    template <class T, class pS>
    void tofile(types::ndarray<T, pS> const &expr, types::str const &file_name,
                types::str const &sep, types::str const &format)
    {
      if (sep.size() != 0)
        throw types::NotImplementedError("Sep input is not implemented yet, should be left empty");
      if (format.size() != 0)
        throw types::NotImplementedError(
            "Format input is not implemented yet, should be left empty");
      std::ofstream fs;
      fs.open(file_name.c_str(), std::ofstream::out | std::ofstream::binary);
      if (fs.rdstate() != std::ofstream::goodbit) {
        throw types::FileNotFoundError("Could not open file " + file_name);
      }
      fs.write((char *)expr.buffer, sizeof(T) * expr.flat_size());
    }
    NUMPY_EXPR_TO_NDARRAY0_IMPL(tofile);
  } // namespace ndarray
} // namespace numpy

PYTHONIC_NS_END

#endif
