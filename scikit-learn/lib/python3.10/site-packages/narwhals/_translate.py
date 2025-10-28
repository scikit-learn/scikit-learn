"""[Protocols] defining conversion methods between representations.

These come in 3 flavors and are [generic] to promote reuse.

The following examples use the placeholder types `Narwhal` and `Other`:
- `Narwhal`: some class written in `narwhals`.
- `Other`: any other class, could be native, compliant, or a builtin.

## `To<Other>`
When we want to convert or unwrap a `Narwhal` into an `Other`,
we provide an **instance** method:

    ToOtherT_co = TypeVar("ToOtherT_co", covariant=True)

    class ToOther(Protocol[ToOtherT_co]):
        def to_other(self, *args: Any, **kwds: Any) -> ToOtherT_co: ...

- `*args`, `**kwds` are defined to be *permissive* and allow a wider set of signatures when implementing.
  - In most cases, they are unused.
  - But come in handy when adapting an [upstream signature].
- We use a  **covariant** `TypeVar`.

## `From<Other>`
But what if we have `Other` and want to do the reverse?

Our `Narwhal` will need to provide a `@classmethod`:

    FromOtherT_contra = TypeVar("FromOtherT_contra", contravariant=True)

    class FromOther(Protocol[FromOtherT_contra]):
        @classmethod
        def from_other(cls, data: FromOtherT_contra, *args: Any, **kwds: Any) -> Self: ...

- `*args`, `**kwds` serve a similar purpose as before, but are much more frequently used.
- We've added a **required** [positional-only] parameter `data` which will always be passed `Other`.
  - This removes the name from the contract of the protocol.
  - Implementations are free to use something more descriptive for documentation purposes.
- We use a  **contravariant** `TypeVar`.

## `<Other>Convertible`
Combining our `to_` and `from_` methods allows us to convert in both directions `Narwhal` <-> `Other`:

    class OtherConvertible(
        ToOther[ToOtherT_co],
        FromOther[FromOtherT_contra],
        Protocol[ToOtherT_co, FromOtherT_contra],
    ): ...

## See Also
Variance of `TypeVar`(s) can be tricky to wrap your head around.

To learn more see [moist], [dry], or [even drier] - depending on how deep you wanna go.

[Protocols]: https://typing.python.org/en/latest/spec/protocol.html
[generic]: https://typing.python.org/en/latest/spec/generics.html
[upstream signature]: https://numpy.org/doc/stable/user/basics.interoperability.html#the-array-method
[positional-only]: https://peps.python.org/pep-0570/
[moist]: https://mypy.readthedocs.io/en/stable/generics.html#variance-of-generic-types
[dry]: https://typing.python.org/en/latest/spec/generics.html#variance
[even drier]: https://en.wikipedia.org/wiki/Covariance_and_contravariance_%28computer_science%29
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import TYPE_CHECKING, Any, Protocol

from narwhals._typing_compat import TypeVar

if TYPE_CHECKING:
    import pyarrow as pa
    from typing_extensions import Self, TypeAlias, TypeIs


class ArrowStreamExportable(Protocol):
    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object: ...


ToNumpyT_co = TypeVar("ToNumpyT_co", covariant=True)
FromNumpyDT_contra = TypeVar(
    "FromNumpyDT_contra", contravariant=True, default=ToNumpyT_co
)
FromNumpyT_contra = TypeVar("FromNumpyT_contra", contravariant=True)


class ToNumpy(Protocol[ToNumpyT_co]):
    def to_numpy(self, *args: Any, **kwds: Any) -> ToNumpyT_co: ...


class FromNumpy(Protocol[FromNumpyT_contra]):
    @classmethod
    def from_numpy(cls, data: FromNumpyT_contra, *args: Any, **kwds: Any) -> Self: ...


class NumpyConvertible(
    ToNumpy[ToNumpyT_co],
    FromNumpy[FromNumpyDT_contra],
    Protocol[ToNumpyT_co, FromNumpyDT_contra],
):
    def to_numpy(self, dtype: Any, *, copy: bool | None) -> ToNumpyT_co: ...


FromIterableT_contra = TypeVar("FromIterableT_contra", contravariant=True, default=Any)


class FromIterable(Protocol[FromIterableT_contra]):
    @classmethod
    def from_iterable(
        cls, data: Iterable[FromIterableT_contra], *args: Any, **kwds: Any
    ) -> Self: ...


ToDictDT_co = TypeVar(
    "ToDictDT_co", bound=Mapping[str, Any], covariant=True, default="dict[str, Any]"
)
FromDictDT_contra = TypeVar(
    "FromDictDT_contra",
    bound=Mapping[str, Any],
    contravariant=True,
    default=Mapping[str, Any],
)


class ToDict(Protocol[ToDictDT_co]):
    def to_dict(self, *args: Any, **kwds: Any) -> ToDictDT_co: ...


class FromDict(Protocol[FromDictDT_contra]):
    @classmethod
    def from_dict(cls, data: FromDictDT_contra, *args: Any, **kwds: Any) -> Self: ...


class DictConvertible(
    ToDict[ToDictDT_co],
    FromDict[FromDictDT_contra],
    Protocol[ToDictDT_co, FromDictDT_contra],
): ...


IntoArrowTable: TypeAlias = "ArrowStreamExportable | pa.Table"
"""An object supporting the [Arrow PyCapsule Interface], or a native [`pyarrow.Table`].

[Arrow PyCapsule Interface]: https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html#arrowstream-export
[`pyarrow.Table`]: https://arrow.apache.org/docs/python/generated/pyarrow.Table.html
"""
ToArrowT_co = TypeVar("ToArrowT_co", covariant=True)
FromArrowDT_contra = TypeVar(
    "FromArrowDT_contra", contravariant=True, default=IntoArrowTable
)


class ToArrow(Protocol[ToArrowT_co]):
    def to_arrow(self, *args: Any, **kwds: Any) -> ToArrowT_co: ...


class FromArrow(Protocol[FromArrowDT_contra]):
    @classmethod
    def from_arrow(cls, data: FromArrowDT_contra, *args: Any, **kwds: Any) -> Self: ...


class ArrowConvertible(
    ToArrow[ToArrowT_co],
    FromArrow[FromArrowDT_contra],
    Protocol[ToArrowT_co, FromArrowDT_contra],
): ...


FromNativeT = TypeVar("FromNativeT")


class FromNative(Protocol[FromNativeT]):
    @classmethod
    def from_native(cls, data: FromNativeT, *args: Any, **kwds: Any) -> Self: ...
    @staticmethod
    def _is_native(obj: FromNativeT | Any, /) -> TypeIs[FromNativeT]:
        """Return `True` if `obj` can be passed to `from_native`."""
        ...


ToNarwhalsT_co = TypeVar("ToNarwhalsT_co", covariant=True)


class ToNarwhals(Protocol[ToNarwhalsT_co]):
    def to_narwhals(self) -> ToNarwhalsT_co:
        """Convert into public representation."""
        ...
