from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal

from narwhals.dependencies import is_narwhals_series

if TYPE_CHECKING:
    from typing_extensions import Never, TypeAlias

# NOTE: This alias is created to facilitate autocomplete. Feel free to extend it as
# you please when adding a new feature.
# See: https://github.com/narwhals-dev/narwhals/pull/2983#discussion_r2337548736
SeriesDetail: TypeAlias = Literal[
    "implementation mismatch",
    "length mismatch",
    "dtype mismatch",
    "name mismatch",
    "null value mismatch",
    "exact value mismatch",
    "values not within tolerance",
    "nested value mismatch",
]


def raise_assertion_error(
    objects: str, detail: str, left: Any, right: Any, *, cause: Exception | None = None
) -> Never:
    """Raise a detailed assertion error."""
    __tracebackhide__ = True

    trailing_left = "\n" if is_narwhals_series(left) else " "
    trailing_right = "\n" if is_narwhals_series(right) else " "

    msg = (
        f"{objects} are different ({detail})\n"
        f"[left]:{trailing_left}{left}\n"
        f"[right]:{trailing_right}{right}"
    )
    raise AssertionError(msg) from cause


def raise_series_assertion_error(
    detail: SeriesDetail, left: Any, right: Any, *, cause: Exception | None = None
) -> Never:
    raise_assertion_error("Series", detail, left, right, cause=cause)
