#ifndef PYTHONIC_NUMPY_FROMFILE_HPP
#define PYTHONIC_NUMPY_FROMFILE_HPP

#include "pythonic/include/numpy/fromfile.hpp"

#include "pythonic/builtins/FileNotFoundError.hpp"
#include "pythonic/builtins/NotImplementedError.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"
#include <fstream>

#include <limits>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype>
  types::ndarray<typename dtype::type, types::pshape<long>>
  fromfile(types::str const &file_name, dtype d, long count, types::str const &sep, long offset)
  {
    if (sep.size() != 0)
      throw types::NotImplementedError("Sep input is not implemented yet, should be left empty");
    std::fstream fs;
    fs.open(file_name.c_str(), std::fstream::in | std::fstream::binary);
    if (fs.rdstate() != std::fstream::goodbit) {
      throw types::FileNotFoundError("Could not find file " + file_name);
    }
    fs.seekp(offset, std::fstream::beg);
    auto n1 = fs.tellp();
    fs.seekp(0, std::fstream::end);
    long maxCount = (fs.tellp() - n1) / sizeof(typename dtype::type);
    fs.seekp(offset, std::fstream::beg);
    if (count < 0) {
      count = maxCount;
    } else if (count > maxCount) {
      count = maxCount;
    }

    types::ndarray<typename dtype::type, types::pshape<long>> res(types::pshape<long>{count},
                                                                  types::none_type{});
    fs.read((char *)res.buffer, sizeof(typename dtype::type) * count);
    return res;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
