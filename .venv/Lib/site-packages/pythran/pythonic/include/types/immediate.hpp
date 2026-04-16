#ifndef PYTHONIC_INCLUDE_TYPES_IMMEDIATE_HPP
#define PYTHONIC_INCLUDE_TYPES_IMMEDIATE_HPP

PYTHONIC_NS_BEGIN

namespace types
{

  template <class T, T Val>
  struct immediate {
    immediate() = default;
    immediate(immediate const &) = default;
    immediate(immediate &&) = default;

    operator T() const
    {
      return Val;
    }

    template <class U, U Wal, class _ = std::enable_if_t<Val == (T)Wal, void>>
    immediate(std::integral_constant<U, Wal>)
    {
    }
  };

  using true_immediate = immediate<bool, true>;
  using false_immediate = immediate<bool, false>;
} // namespace types

PYTHONIC_NS_END

#endif
