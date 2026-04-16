#ifndef PYTHONIC_BUILTIN_XRANGE_HPP
#define PYTHONIC_BUILTIN_XRANGE_HPP

#include "pythonic/include/builtins/xrange.hpp"

#include "pythonic/utils/functor.hpp"

#include <iterator>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace
  {

    long _init_last(long _begin, long _end, long _step)
    {
      if (_step > 0)
        return _begin + std::max(0L, _step * ((_end - _begin + _step - 1) / _step));
      else
        return _begin + std::min(0L, _step * ((_end - _begin + _step + 1) / _step));
    }
  } // namespace

  xrange_iterator::xrange_iterator(long v, long s) : value_(v), step_(s)
  {
  }

  long xrange_iterator::operator*() const
  {
    return value_;
  }

  xrange_iterator &xrange_iterator::operator++()
  {
    value_ += step_;
    return *this;
  }

  xrange_iterator xrange_iterator::operator++(int)
  {
    xrange_iterator self(*this);
    value_ += step_;
    return self;
  }

  xrange_iterator &xrange_iterator::operator+=(long n)
  {
    value_ += step_ * n;
    return *this;
  }

  bool xrange_iterator::operator!=(xrange_iterator const &other) const
  {
    return value_ != other.value_;
  }

  bool xrange_iterator::operator==(xrange_iterator const &other) const
  {
    return value_ == other.value_;
  }

  bool xrange_iterator::operator<(xrange_iterator const &other) const
  {
    return step_ * value_ < step_ * other.value_;
  }

  long xrange_iterator::operator-(xrange_iterator const &other) const
  {
    return (value_ - other.value_) / step_;
  }

  xrange::xrange(long b, long e, long s) : begin_(b), end_(_init_last(b, e, s)), step_(s)
  {
  }

  xrange::xrange(long e) : begin_(0), end_(e), step_(1)
  {
  }

  xrange_iterator xrange::begin() const
  {
    return xrange_iterator(begin_, step_);
  }

  xrange_iterator xrange::end() const
  {
    return xrange_iterator(end_, step_);
  }

  typename xrange::reverse_iterator xrange::rbegin() const
  {
    return {end_ - step_, -step_};
  }

  typename xrange::reverse_iterator xrange::rend() const
  {
    return {begin_ - step_, -step_};
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
