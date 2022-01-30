"""Data structures configuring Black behavior.

Mostly around Python language feature support per version and Black configuration
chosen by the user.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, Set

from black.const import DEFAULT_LINE_LENGTH


class TargetVersion(Enum):
    PY27 = 2
    PY33 = 3
    PY34 = 4
    PY35 = 5
    PY36 = 6
    PY37 = 7
    PY38 = 8
    PY39 = 9

    def is_python2(self) -> bool:
        return self is TargetVersion.PY27


class Feature(Enum):
    # All string literals are unicode
    UNICODE_LITERALS = 1
    F_STRINGS = 2
    NUMERIC_UNDERSCORES = 3
    TRAILING_COMMA_IN_CALL = 4
    TRAILING_COMMA_IN_DEF = 5
    # The following two feature-flags are mutually exclusive, and exactly one should be
    # set for every version of python.
    ASYNC_IDENTIFIERS = 6
    ASYNC_KEYWORDS = 7
    ASSIGNMENT_EXPRESSIONS = 8
    POS_ONLY_ARGUMENTS = 9
    RELAXED_DECORATORS = 10
    FORCE_OPTIONAL_PARENTHESES = 50


VERSION_TO_FEATURES: Dict[TargetVersion, Set[Feature]] = {
    TargetVersion.PY27: {Feature.ASYNC_IDENTIFIERS},
    TargetVersion.PY33: {Feature.UNICODE_LITERALS, Feature.ASYNC_IDENTIFIERS},
    TargetVersion.PY34: {Feature.UNICODE_LITERALS, Feature.ASYNC_IDENTIFIERS},
    TargetVersion.PY35: {
        Feature.UNICODE_LITERALS,
        Feature.TRAILING_COMMA_IN_CALL,
        Feature.ASYNC_IDENTIFIERS,
    },
    TargetVersion.PY36: {
        Feature.UNICODE_LITERALS,
        Feature.F_STRINGS,
        Feature.NUMERIC_UNDERSCORES,
        Feature.TRAILING_COMMA_IN_CALL,
        Feature.TRAILING_COMMA_IN_DEF,
        Feature.ASYNC_IDENTIFIERS,
    },
    TargetVersion.PY37: {
        Feature.UNICODE_LITERALS,
        Feature.F_STRINGS,
        Feature.NUMERIC_UNDERSCORES,
        Feature.TRAILING_COMMA_IN_CALL,
        Feature.TRAILING_COMMA_IN_DEF,
        Feature.ASYNC_KEYWORDS,
    },
    TargetVersion.PY38: {
        Feature.UNICODE_LITERALS,
        Feature.F_STRINGS,
        Feature.NUMERIC_UNDERSCORES,
        Feature.TRAILING_COMMA_IN_CALL,
        Feature.TRAILING_COMMA_IN_DEF,
        Feature.ASYNC_KEYWORDS,
        Feature.ASSIGNMENT_EXPRESSIONS,
        Feature.POS_ONLY_ARGUMENTS,
    },
    TargetVersion.PY39: {
        Feature.UNICODE_LITERALS,
        Feature.F_STRINGS,
        Feature.NUMERIC_UNDERSCORES,
        Feature.TRAILING_COMMA_IN_CALL,
        Feature.TRAILING_COMMA_IN_DEF,
        Feature.ASYNC_KEYWORDS,
        Feature.ASSIGNMENT_EXPRESSIONS,
        Feature.RELAXED_DECORATORS,
        Feature.POS_ONLY_ARGUMENTS,
    },
}


def supports_feature(target_versions: Set[TargetVersion], feature: Feature) -> bool:
    return all(feature in VERSION_TO_FEATURES[version] for version in target_versions)


@dataclass
class Mode:
    target_versions: Set[TargetVersion] = field(default_factory=set)
    line_length: int = DEFAULT_LINE_LENGTH
    string_normalization: bool = True
    is_pyi: bool = False
    magic_trailing_comma: bool = True
    experimental_string_processing: bool = False

    def get_cache_key(self) -> str:
        if self.target_versions:
            version_str = ",".join(
                str(version.value)
                for version in sorted(self.target_versions, key=lambda v: v.value)
            )
        else:
            version_str = "-"
        parts = [
            version_str,
            str(self.line_length),
            str(int(self.string_normalization)),
            str(int(self.is_pyi)),
            str(int(self.magic_trailing_comma)),
            str(int(self.experimental_string_processing)),
        ]
        return ".".join(parts)
