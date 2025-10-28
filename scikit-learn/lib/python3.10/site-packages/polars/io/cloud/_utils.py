from __future__ import annotations

from pathlib import Path
from typing import Any, Generic, TypeVar

from polars._typing import PartitioningScheme
from polars._utils.various import is_path_or_str_sequence

T = TypeVar("T")


class NoPickleOption(Generic[T]):
    """
    Wrapper that does not pickle the wrapped value.

    This wrapper will unpickle to contain a None. Used for cached values.
    """

    def __init__(self, opt_value: T | None = None) -> None:
        self._opt_value = opt_value

    def get(self) -> T | None:
        return self._opt_value

    def set(self, value: T | None) -> None:
        self._opt_value = value

    def __getstate__(self) -> tuple[()]:
        # Needs to return not-None for `__setstate__()` to be called
        return ()

    def __setstate__(self, _state: tuple[()]) -> None:
        NoPickleOption.__init__(self)


def _first_scan_path(
    source: Any,
) -> str | Path | None:
    if isinstance(source, (str, Path)):
        return source
    elif is_path_or_str_sequence(source) and source:
        return source[0]
    elif isinstance(source, PartitioningScheme):
        return source._base_path

    return None


def _get_path_scheme(path: str | Path) -> str | None:
    path_str = str(path)
    i = path_str.find("://")

    return path_str[:i] if i >= 0 else None


def _is_aws_cloud(*, scheme: str, first_scan_path: str) -> bool:
    if any(scheme == x for x in ["s3", "s3a"]):
        return True

    if scheme == "http" or scheme == "https":
        bucket_end = first_scan_path.find(".s3.")
        region_end = first_scan_path.find(".amazonaws.com/", bucket_end + 4)

        if (
            first_scan_path.find("/", len(scheme) + 3, region_end) > 0
            or "?" in first_scan_path
        ):
            return False

        return 0 < bucket_end < region_end

    return False


def _is_azure_cloud(scheme: str) -> bool:
    return any(scheme == x for x in ["az", "azure", "adl", "abfs", "abfss"])


def _is_gcp_cloud(scheme: str) -> bool:
    return any(scheme == x for x in ["gs", "gcp", "gcs"])
