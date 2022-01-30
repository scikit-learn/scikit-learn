"""Classification of possible errors mypy can detect.

These can be used for filtering specific errors.
"""

from typing import Dict, List
from typing_extensions import Final


# All created error codes are implicitly stored in this list.
all_error_codes: List["ErrorCode"] = []

error_codes: Dict[str, "ErrorCode"] = {}


class ErrorCode:
    def __init__(self, code: str,
                 description: str,
                 category: str,
                 default_enabled: bool = True) -> None:
        self.code = code
        self.description = description
        self.category = category
        self.default_enabled = default_enabled
        error_codes[code] = self

    def __str__(self) -> str:
        return '<ErrorCode {}>'.format(self.code)


ATTR_DEFINED: Final = ErrorCode("attr-defined", "Check that attribute exists", "General")
NAME_DEFINED: Final = ErrorCode("name-defined", "Check that name is defined", "General")
CALL_ARG: Final[ErrorCode] = ErrorCode(
    "call-arg", "Check number, names and kinds of arguments in calls", "General"
)
ARG_TYPE: Final = ErrorCode("arg-type", "Check argument types in calls", "General")
CALL_OVERLOAD: Final = ErrorCode(
    "call-overload", "Check that an overload variant matches arguments", "General"
)
VALID_TYPE: Final = ErrorCode("valid-type", "Check that type (annotation) is valid", "General")
VAR_ANNOTATED: Final = ErrorCode(
    "var-annotated", "Require variable annotation if type can't be inferred", "General"
)
OVERRIDE: Final = ErrorCode(
    "override", "Check that method override is compatible with base class", "General"
)
RETURN: Final[ErrorCode] = ErrorCode(
    "return", "Check that function always returns a value", "General"
)
RETURN_VALUE: Final[ErrorCode] = ErrorCode(
    "return-value", "Check that return value is compatible with signature", "General"
)
ASSIGNMENT: Final = ErrorCode(
    "assignment", "Check that assigned value is compatible with target", "General"
)
TYPE_ARG: Final = ErrorCode("type-arg", "Check that generic type arguments are present", "General")
TYPE_VAR: Final = ErrorCode("type-var", "Check that type variable values are valid", "General")
UNION_ATTR: Final = ErrorCode(
    "union-attr", "Check that attribute exists in each item of a union", "General"
)
INDEX: Final = ErrorCode("index", "Check indexing operations", "General")
OPERATOR: Final = ErrorCode("operator", "Check that operator is valid for operands", "General")
LIST_ITEM: Final = ErrorCode(
    "list-item", "Check list items in a list expression [item, ...]", "General"
)
DICT_ITEM: Final = ErrorCode(
    "dict-item", "Check dict items in a dict expression {key: value, ...}", "General"
)
TYPEDDICT_ITEM: Final = ErrorCode(
    "typeddict-item", "Check items when constructing TypedDict", "General"
)
HAS_TYPE: Final = ErrorCode(
    "has-type", "Check that type of reference can be determined", "General"
)
IMPORT: Final = ErrorCode(
    "import", "Require that imported module can be found or has stubs", "General"
)
NO_REDEF: Final = ErrorCode("no-redef", "Check that each name is defined once", "General")
FUNC_RETURNS_VALUE: Final = ErrorCode(
    "func-returns-value", "Check that called function returns a value in value context", "General"
)
ABSTRACT: Final = ErrorCode(
    "abstract", "Prevent instantiation of classes with abstract attributes", "General"
)
VALID_NEWTYPE: Final = ErrorCode(
    "valid-newtype", "Check that argument 2 to NewType is valid", "General"
)
STRING_FORMATTING: Final = ErrorCode(
    "str-format", "Check that string formatting/interpolation is type-safe", "General"
)
STR_BYTES_PY3: Final = ErrorCode(
    "str-bytes-safe", "Warn about dangerous coercions related to bytes and string types", "General"
)
EXIT_RETURN: Final = ErrorCode(
    "exit-return", "Warn about too general return type for '__exit__'", "General"
)

# These error codes aren't enabled by default.
NO_UNTYPED_DEF: Final[ErrorCode] = ErrorCode(
    "no-untyped-def", "Check that every function has an annotation", "General"
)
NO_UNTYPED_CALL: Final = ErrorCode(
    "no-untyped-call",
    "Disallow calling functions without type annotations from annotated functions",
    "General",
)
REDUNDANT_CAST: Final = ErrorCode(
    "redundant-cast", "Check that cast changes type of expression", "General"
)
COMPARISON_OVERLAP: Final = ErrorCode(
    "comparison-overlap", "Check that types in comparisons and 'in' expressions overlap", "General"
)
NO_ANY_UNIMPORTED: Final = ErrorCode(
    "no-any-unimported", 'Reject "Any" types from unfollowed imports', "General"
)
NO_ANY_RETURN: Final = ErrorCode(
    "no-any-return",
    'Reject returning value with "Any" type if return type is not "Any"',
    "General",
)
UNREACHABLE: Final = ErrorCode(
    "unreachable", "Warn about unreachable statements or expressions", "General"
)
REDUNDANT_EXPR: Final = ErrorCode(
    "redundant-expr", "Warn about redundant expressions", "General", default_enabled=False
)
TRUTHY_BOOL: Final[ErrorCode] = ErrorCode(
    "truthy-bool",
    "Warn about expressions that could always evaluate to true in boolean contexts",
    "General",
    default_enabled=False,
)
NAME_MATCH: Final = ErrorCode(
    "name-match", "Check that type definition has consistent naming", "General"
)


# Syntax errors are often blocking.
SYNTAX: Final = ErrorCode("syntax", "Report syntax errors", "General")

# This is a catch-all for remaining uncategorized errors.
MISC: Final = ErrorCode("misc", "Miscellaneous other checks", "General")
