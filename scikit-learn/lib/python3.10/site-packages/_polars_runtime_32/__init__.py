# ruff: noqa: TID252
try:
    from .build_feature_flags import BUILD_FEATURE_FLAGS
except ImportError:
    BUILD_FEATURE_FLAGS = ""

__all__ = ["BUILD_FEATURE_FLAGS"]
