#ifndef PYTHONIC_NUMPY_CORRELATE_HPP
#define PYTHONIC_NUMPY_CORRELATE_HPP

#include "pythonic/include/numpy/correlate.hpp"
#include "pythonic/numpy/asarray.hpp"
#include "pythonic/numpy/conjugate.hpp"
#include "pythonic/numpy/dot.hpp"
#include "pythonic/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class A, class B>
  types::ndarray<typename A::dtype, types::pshape<long>>
  do_correlate(A const &inA, B const &inB, types::str const &type, int out_inc)
  // out_inc is used to indicate the inputs were swapped, which means that the
  // output must be time reversed and conjugated
  {
    auto shapeA = sutils::getshape(inA);
    auto shapeB = sutils::getshape(inB);

    long NA = shapeA[0];
    long NB = shapeB[0];

    using out_type = typename __combined<typename A::dtype, typename B::dtype>::type;
    // At this point, handling views would slow things down tremendously
    auto inA_ = functor::asarray{}(inA);
    auto inB_ = functor::asarray{}(inB);

    auto outN = 0;
    int iLeft;
    if (type == "full") {
      outN = NA + NB - 1;
      iLeft = -NB + 1;
    } else if (type == "valid") {
      outN = NA - NB + 1;
      iLeft = 0;
    } else {
      assert(type == "same" && "valid type");
      outN = NA;
      iLeft = -NB + 1 + (NB - 1) / 2;
    }
    // We need outN output values, no matter what.
    int iRight = iLeft + outN;

    // Allocate output array
    types::ndarray<out_type, types::pshape<long>> out = {outN, out_type()};
    out_type *out_ptr = (out_type *)out.buffer;
    // if out_inc is -1, we reverse the output.
    if (out_inc == -1)
      out_ptr += outN - 1;

    // For small correlations, numpy uses small_correlate, far more efficient.
    // see numpy/core/src/multiarray/arraytypes.c.src

    if (out_inc == 1) {
      // Incomplete overlap left
      for (int i = iLeft; i < 0; i++, out_ptr++) {
        *out_ptr = numpy::dot(inA_(types::fast_contiguous_slice(0, NB + i)),
                              inB_(types::fast_contiguous_slice(-i, NB)));
      }
      // Complete overlap middle
      for (int i = 0; i <= NA - NB; i++, out_ptr++) {
        *out_ptr = numpy::dot(inA_(types::fast_contiguous_slice(i, i + NB)),
                              inB_(types::fast_contiguous_slice(0, NB)));
      }
      // Incomplete overlap right.
      for (int i = NA - NB + 1; i < iRight; i++, out_ptr++) {
        *out_ptr = numpy::dot(inA_(types::fast_contiguous_slice(i, NA)),
                              inB_(types::fast_contiguous_slice(0, NA - i)));
      }
    } else {
      // Incomplete overlap left
      for (int i = iLeft; i < 0; i++, out_ptr += out_inc) {
        *out_ptr = wrapper::conjugate(numpy::dot(inA_(types::fast_contiguous_slice(0, NB + i)),
                                                 inB_(types::fast_contiguous_slice(-i, NB))));
      }
      // Complete overlap middle
      for (int i = 0; i <= NA - NB; i++, out_ptr += out_inc) {
        *out_ptr = wrapper::conjugate(numpy::dot(inA_(types::fast_contiguous_slice(i, i + NB)),
                                                 inB_(types::fast_contiguous_slice(0, NB))));
      }
      // Incomplete overlap right.
      for (int i = NA - NB + 1; i < iRight; i++, out_ptr += out_inc) {
        *out_ptr = wrapper::conjugate(numpy::dot(inA_(types::fast_contiguous_slice(i, NA)),
                                                 inB_(types::fast_contiguous_slice(0, NA - i))));
      }
    }

    return out;
  }

  template <class A, class B>
  types::ndarray<typename A::dtype, types::pshape<long>> correlate(A const &inA, B const &inB,
                                                                   types::str const &type)
  {
    long NA = inA.template shape<0>();
    long NB = inB.template shape<0>();
    // If inB is longer than inA, swap them, but time-reverse and conjugate the
    // output (-1 flag)
    if (NA > NB) {
      auto inB_conj = functor::conjugate{}(inB);
      return do_correlate(inA, inB_conj, type, 1);
    } else {
      auto inA_conj = functor::conjugate{}(inA);
      return do_correlate(inB, inA_conj, type, -1);
    }
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(correlate)
} // namespace numpy
PYTHONIC_NS_END

#endif
