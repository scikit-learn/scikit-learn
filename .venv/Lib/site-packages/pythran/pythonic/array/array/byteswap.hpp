#ifndef PYTHONIC_ARRAY_ARRAY_BYTESWAP_HPP
#define PYTHONIC_ARRAY_ARRAY_BYTESWAP_HPP

#include "pythonic/include/array/array/byteswap.hpp"
#include "pythonic/utils/functor.hpp"

#include <byteswap.h>

PYTHONIC_NS_BEGIN

namespace array
{

  namespace array
  {
    inline void byteswap(char *buffer, size_t n, std::integral_constant<unsigned, 2>)
    {
      auto *ibuffer = reinterpret_cast<uint16_t *>(buffer);
      for (size_t i = 0; i < n; i++)
        ibuffer[i] = bswap_16(ibuffer[i]);
    }
    inline void byteswap(char *buffer, size_t n, std::integral_constant<unsigned, 4>)
    {
      auto *ibuffer = reinterpret_cast<uint32_t *>(buffer);
      for (size_t i = 0; i < n; i++)
        ibuffer[i] = bswap_32(ibuffer[i]);
    }
    inline void byteswap(char *buffer, size_t n, std::integral_constant<unsigned, 8>)
    {
      auto *ibuffer = reinterpret_cast<uint64_t *>(buffer);
      for (size_t i = 0; i < n; i++)
        ibuffer[i] = bswap_64(ibuffer[i]);
    }

    template <class T>
    types::none_type byteswap(types::array<T> &seq)
    {
      byteswap(reinterpret_cast<char *>(seq.data()), seq.size(),
               std::integral_constant<unsigned, sizeof(T)>{});
      return {};
    }

  } // namespace array
} // namespace array
PYTHONIC_NS_END
#endif
