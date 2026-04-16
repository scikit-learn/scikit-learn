#ifndef PYTHONIC_INCLUDE_TYPES_ATTR_HPP
#define PYTHONIC_INCLUDE_TYPES_ATTR_HPP

PYTHONIC_NS_BEGIN

namespace types
{

  namespace attr
  {
    /* exception attributes */
    struct ARGS {
    };
    struct ERRNO {
    };
    struct STRERROR {
    };
    struct FILENAME {
    };
    /* complex attributes */
    struct REAL {
    };
    struct IMAG {
    };
    /* file attributes */
    struct CLOSED {
    };
    struct MODE {
    };
    struct NAME {
    };
    struct NEWLINES {
    };
    /* fileinfo attributes */
    struct EPS {
    };
    /* ndarray attributes */
    struct SHAPE {
    };
    struct NDIM {
    };
    struct STRIDES {
    };
    struct SIZE {
    };
    struct ITEMSIZE {
    };
    struct NBYTES {
    };
    struct FLAT {
    };
    struct DTYPE {
    };
    struct T {
    };
    /* slice attributes */
    struct START {
    };
    struct STOP {
    };
    struct STEP {
    };
    /* */
  } // namespace attr
} // namespace types
PYTHONIC_NS_END

#endif
