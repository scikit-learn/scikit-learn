import sys
from typing import (
    IO,
    Any,
    Callable,
    Generator,
    Generic,
    Iterable,
    NoReturn,
    Pattern,
    Protocol,
    Sequence,
    Tuple,
    Type,
    TypeVar,
    overload,
)

_T = TypeVar("_T")
_ActionT = TypeVar("_ActionT", bound=Action)
_ArgumentParserT = TypeVar("_ArgumentParserT", bound=ArgumentParser)
_N = TypeVar("_N")

ONE_OR_MORE: str
OPTIONAL: str
PARSER: str
REMAINDER: str
SUPPRESS: str
ZERO_OR_MORE: str
_UNRECOGNIZED_ARGS_ATTR: str  # undocumented

class ArgumentError(Exception):
    argument_name: str | None
    message: str
    def __init__(self, argument: Action | None, message: str) -> None: ...

# undocumented
class _AttributeHolder:
    def _get_kwargs(self) -> list[tuple[str, Any]]: ...
    def _get_args(self) -> list[Any]: ...

# undocumented
class _ActionsContainer:
    description: str | None
    prefix_chars: str
    argument_default: Any
    conflict_handler: str

    _registries: dict[str, dict[Any, Any]]
    _actions: list[Action]
    _option_string_actions: dict[str, Action]
    _action_groups: list[_ArgumentGroup]
    _mutually_exclusive_groups: list[_MutuallyExclusiveGroup]
    _defaults: dict[str, Any]
    _negative_number_matcher: Pattern[str]
    _has_negative_number_optionals: list[bool]
    def __init__(self, description: str | None, prefix_chars: str, argument_default: Any, conflict_handler: str) -> None: ...
    def register(self, registry_name: str, value: Any, object: Any) -> None: ...
    def _registry_get(self, registry_name: str, value: Any, default: Any = ...) -> Any: ...
    def set_defaults(self, **kwargs: Any) -> None: ...
    def get_default(self, dest: str) -> Any: ...
    def add_argument(
        self,
        *name_or_flags: str,
        action: str | Type[Action] = ...,
        nargs: int | str = ...,
        const: Any = ...,
        default: Any = ...,
        type: Callable[[str], _T] | Callable[[str], _T] | FileType = ...,
        choices: Iterable[_T] | None = ...,
        required: bool = ...,
        help: str | None = ...,
        metavar: str | Tuple[str, ...] | None = ...,
        dest: str | None = ...,
        version: str = ...,
        **kwargs: Any,
    ) -> Action: ...
    def add_argument_group(self, *args: Any, **kwargs: Any) -> _ArgumentGroup: ...
    def add_mutually_exclusive_group(self, **kwargs: Any) -> _MutuallyExclusiveGroup: ...
    def _add_action(self, action: _ActionT) -> _ActionT: ...
    def _remove_action(self, action: Action) -> None: ...
    def _add_container_actions(self, container: _ActionsContainer) -> None: ...
    def _get_positional_kwargs(self, dest: str, **kwargs: Any) -> dict[str, Any]: ...
    def _get_optional_kwargs(self, *args: Any, **kwargs: Any) -> dict[str, Any]: ...
    def _pop_action_class(self, kwargs: Any, default: Type[Action] | None = ...) -> Type[Action]: ...
    def _get_handler(self) -> Callable[[Action, Iterable[tuple[str, Action]]], Any]: ...
    def _check_conflict(self, action: Action) -> None: ...
    def _handle_conflict_error(self, action: Action, conflicting_actions: Iterable[tuple[str, Action]]) -> NoReturn: ...
    def _handle_conflict_resolve(self, action: Action, conflicting_actions: Iterable[tuple[str, Action]]) -> None: ...

class _FormatterClass(Protocol):
    def __call__(self, prog: str) -> HelpFormatter: ...

class ArgumentParser(_AttributeHolder, _ActionsContainer):
    prog: str
    usage: str | None
    epilog: str | None
    formatter_class: _FormatterClass
    fromfile_prefix_chars: str | None
    add_help: bool
    allow_abbrev: bool

    # undocumented
    _positionals: _ArgumentGroup
    _optionals: _ArgumentGroup
    _subparsers: _ArgumentGroup | None

    if sys.version_info >= (3, 9):
        def __init__(
            self,
            prog: str | None = ...,
            usage: str | None = ...,
            description: str | None = ...,
            epilog: str | None = ...,
            parents: Sequence[ArgumentParser] = ...,
            formatter_class: _FormatterClass = ...,
            prefix_chars: str = ...,
            fromfile_prefix_chars: str | None = ...,
            argument_default: Any = ...,
            conflict_handler: str = ...,
            add_help: bool = ...,
            allow_abbrev: bool = ...,
            exit_on_error: bool = ...,
        ) -> None: ...
    else:
        def __init__(
            self,
            prog: str | None = ...,
            usage: str | None = ...,
            description: str | None = ...,
            epilog: str | None = ...,
            parents: Sequence[ArgumentParser] = ...,
            formatter_class: _FormatterClass = ...,
            prefix_chars: str = ...,
            fromfile_prefix_chars: str | None = ...,
            argument_default: Any = ...,
            conflict_handler: str = ...,
            add_help: bool = ...,
            allow_abbrev: bool = ...,
        ) -> None: ...
    # The type-ignores in these overloads should be temporary.  See:
    # https://github.com/python/typeshed/pull/2643#issuecomment-442280277
    @overload
    def parse_args(self, args: Sequence[str] | None = ...) -> Namespace: ...
    @overload
    def parse_args(self, args: Sequence[str] | None, namespace: None) -> Namespace: ...  # type: ignore
    @overload
    def parse_args(self, args: Sequence[str] | None, namespace: _N) -> _N: ...
    @overload
    def parse_args(self, *, namespace: None) -> Namespace: ...  # type: ignore
    @overload
    def parse_args(self, *, namespace: _N) -> _N: ...
    if sys.version_info >= (3, 7):
        @overload
        def add_subparsers(
            self: _ArgumentParserT,
            *,
            title: str = ...,
            description: str | None = ...,
            prog: str = ...,
            action: Type[Action] = ...,
            option_string: str = ...,
            dest: str | None = ...,
            required: bool = ...,
            help: str | None = ...,
            metavar: str | None = ...,
        ) -> _SubParsersAction[_ArgumentParserT]: ...
        @overload
        def add_subparsers(
            self,
            *,
            title: str = ...,
            description: str | None = ...,
            prog: str = ...,
            parser_class: Type[_ArgumentParserT] = ...,
            action: Type[Action] = ...,
            option_string: str = ...,
            dest: str | None = ...,
            required: bool = ...,
            help: str | None = ...,
            metavar: str | None = ...,
        ) -> _SubParsersAction[_ArgumentParserT]: ...
    else:
        @overload
        def add_subparsers(
            self: _ArgumentParserT,
            *,
            title: str = ...,
            description: str | None = ...,
            prog: str = ...,
            action: Type[Action] = ...,
            option_string: str = ...,
            dest: str | None = ...,
            help: str | None = ...,
            metavar: str | None = ...,
        ) -> _SubParsersAction[_ArgumentParserT]: ...
        @overload
        def add_subparsers(
            self,
            *,
            title: str = ...,
            description: str | None = ...,
            prog: str = ...,
            parser_class: Type[_ArgumentParserT] = ...,
            action: Type[Action] = ...,
            option_string: str = ...,
            dest: str | None = ...,
            help: str | None = ...,
            metavar: str | None = ...,
        ) -> _SubParsersAction[_ArgumentParserT]: ...
    def print_usage(self, file: IO[str] | None = ...) -> None: ...
    def print_help(self, file: IO[str] | None = ...) -> None: ...
    def format_usage(self) -> str: ...
    def format_help(self) -> str: ...
    def parse_known_args(
        self, args: Sequence[str] | None = ..., namespace: Namespace | None = ...
    ) -> tuple[Namespace, list[str]]: ...
    def convert_arg_line_to_args(self, arg_line: str) -> list[str]: ...
    def exit(self, status: int = ..., message: str | None = ...) -> NoReturn: ...
    def error(self, message: str) -> NoReturn: ...
    if sys.version_info >= (3, 7):
        def parse_intermixed_args(self, args: Sequence[str] | None = ..., namespace: Namespace | None = ...) -> Namespace: ...
        def parse_known_intermixed_args(
            self, args: Sequence[str] | None = ..., namespace: Namespace | None = ...
        ) -> tuple[Namespace, list[str]]: ...
    # undocumented
    def _get_optional_actions(self) -> list[Action]: ...
    def _get_positional_actions(self) -> list[Action]: ...
    def _parse_known_args(self, arg_strings: list[str], namespace: Namespace) -> tuple[Namespace, list[str]]: ...
    def _read_args_from_files(self, arg_strings: list[str]) -> list[str]: ...
    def _match_argument(self, action: Action, arg_strings_pattern: str) -> int: ...
    def _match_arguments_partial(self, actions: Sequence[Action], arg_strings_pattern: str) -> list[int]: ...
    def _parse_optional(self, arg_string: str) -> tuple[Action | None, str, str | None] | None: ...
    def _get_option_tuples(self, option_string: str) -> list[tuple[Action, str, str | None]]: ...
    def _get_nargs_pattern(self, action: Action) -> str: ...
    def _get_values(self, action: Action, arg_strings: list[str]) -> Any: ...
    def _get_value(self, action: Action, arg_string: str) -> Any: ...
    def _check_value(self, action: Action, value: Any) -> None: ...
    def _get_formatter(self) -> HelpFormatter: ...
    def _print_message(self, message: str, file: IO[str] | None = ...) -> None: ...

class HelpFormatter:
    # undocumented
    _prog: str
    _indent_increment: int
    _max_help_position: int
    _width: int
    _current_indent: int
    _level: int
    _action_max_length: int
    _root_section: Any
    _current_section: Any
    _whitespace_matcher: Pattern[str]
    _long_break_matcher: Pattern[str]
    _Section: Type[Any]  # Nested class
    def __init__(self, prog: str, indent_increment: int = ..., max_help_position: int = ..., width: int | None = ...) -> None: ...
    def _indent(self) -> None: ...
    def _dedent(self) -> None: ...
    def _add_item(self, func: Callable[..., str], args: Iterable[Any]) -> None: ...
    def start_section(self, heading: str | None) -> None: ...
    def end_section(self) -> None: ...
    def add_text(self, text: str | None) -> None: ...
    def add_usage(
        self, usage: str | None, actions: Iterable[Action], groups: Iterable[_ArgumentGroup], prefix: str | None = ...
    ) -> None: ...
    def add_argument(self, action: Action) -> None: ...
    def add_arguments(self, actions: Iterable[Action]) -> None: ...
    def format_help(self) -> str: ...
    def _join_parts(self, part_strings: Iterable[str]) -> str: ...
    def _format_usage(
        self, usage: str, actions: Iterable[Action], groups: Iterable[_ArgumentGroup], prefix: str | None
    ) -> str: ...
    def _format_actions_usage(self, actions: Iterable[Action], groups: Iterable[_ArgumentGroup]) -> str: ...
    def _format_text(self, text: str) -> str: ...
    def _format_action(self, action: Action) -> str: ...
    def _format_action_invocation(self, action: Action) -> str: ...
    def _metavar_formatter(self, action: Action, default_metavar: str) -> Callable[[int], Tuple[str, ...]]: ...
    def _format_args(self, action: Action, default_metavar: str) -> str: ...
    def _expand_help(self, action: Action) -> str: ...
    def _iter_indented_subactions(self, action: Action) -> Generator[Action, None, None]: ...
    def _split_lines(self, text: str, width: int) -> list[str]: ...
    def _fill_text(self, text: str, width: int, indent: str) -> str: ...
    def _get_help_string(self, action: Action) -> str | None: ...
    def _get_default_metavar_for_optional(self, action: Action) -> str: ...
    def _get_default_metavar_for_positional(self, action: Action) -> str: ...

class RawDescriptionHelpFormatter(HelpFormatter): ...
class RawTextHelpFormatter(RawDescriptionHelpFormatter): ...
class ArgumentDefaultsHelpFormatter(HelpFormatter): ...
class MetavarTypeHelpFormatter(HelpFormatter): ...

class Action(_AttributeHolder):
    option_strings: Sequence[str]
    dest: str
    nargs: int | str | None
    const: Any
    default: Any
    type: Callable[[str], Any] | FileType | None
    choices: Iterable[Any] | None
    required: bool
    help: str | None
    metavar: str | Tuple[str, ...] | None
    def __init__(
        self,
        option_strings: Sequence[str],
        dest: str,
        nargs: int | str | None = ...,
        const: _T | None = ...,
        default: _T | str | None = ...,
        type: Callable[[str], _T] | Callable[[str], _T] | FileType | None = ...,
        choices: Iterable[_T] | None = ...,
        required: bool = ...,
        help: str | None = ...,
        metavar: str | Tuple[str, ...] | None = ...,
    ) -> None: ...
    def __call__(
        self, parser: ArgumentParser, namespace: Namespace, values: str | Sequence[Any] | None, option_string: str | None = ...
    ) -> None: ...
    if sys.version_info >= (3, 9):
        def format_usage(self) -> str: ...

if sys.version_info >= (3, 9):
    class BooleanOptionalAction(Action):
        def __init__(
            self,
            option_strings: Sequence[str],
            dest: str,
            default: _T | str | None = ...,
            type: Callable[[str], _T] | Callable[[str], _T] | FileType | None = ...,
            choices: Iterable[_T] | None = ...,
            required: bool = ...,
            help: str | None = ...,
            metavar: str | Tuple[str, ...] | None = ...,
        ) -> None: ...

class Namespace(_AttributeHolder):
    def __init__(self, **kwargs: Any) -> None: ...
    def __getattr__(self, name: str) -> Any: ...
    def __setattr__(self, name: str, value: Any) -> None: ...
    def __contains__(self, key: str) -> bool: ...

class FileType:
    # undocumented
    _mode: str
    _bufsize: int
    _encoding: str | None
    _errors: str | None
    def __init__(self, mode: str = ..., bufsize: int = ..., encoding: str | None = ..., errors: str | None = ...) -> None: ...
    def __call__(self, string: str) -> IO[Any]: ...

# undocumented
class _ArgumentGroup(_ActionsContainer):
    title: str | None
    _group_actions: list[Action]
    def __init__(
        self, container: _ActionsContainer, title: str | None = ..., description: str | None = ..., **kwargs: Any
    ) -> None: ...

# undocumented
class _MutuallyExclusiveGroup(_ArgumentGroup):
    required: bool
    _container: _ActionsContainer
    def __init__(self, container: _ActionsContainer, required: bool = ...) -> None: ...

# undocumented
class _StoreAction(Action): ...

# undocumented
class _StoreConstAction(Action):
    def __init__(
        self,
        option_strings: Sequence[str],
        dest: str,
        const: Any,
        default: Any = ...,
        required: bool = ...,
        help: str | None = ...,
        metavar: str | Tuple[str, ...] | None = ...,
    ) -> None: ...

# undocumented
class _StoreTrueAction(_StoreConstAction):
    def __init__(
        self, option_strings: Sequence[str], dest: str, default: bool = ..., required: bool = ..., help: str | None = ...
    ) -> None: ...

# undocumented
class _StoreFalseAction(_StoreConstAction):
    def __init__(
        self, option_strings: Sequence[str], dest: str, default: bool = ..., required: bool = ..., help: str | None = ...
    ) -> None: ...

# undocumented
class _AppendAction(Action): ...

# undocumented
class _AppendConstAction(Action):
    def __init__(
        self,
        option_strings: Sequence[str],
        dest: str,
        const: Any,
        default: Any = ...,
        required: bool = ...,
        help: str | None = ...,
        metavar: str | Tuple[str, ...] | None = ...,
    ) -> None: ...

# undocumented
class _CountAction(Action):
    def __init__(
        self, option_strings: Sequence[str], dest: str, default: Any = ..., required: bool = ..., help: str | None = ...
    ) -> None: ...

# undocumented
class _HelpAction(Action):
    def __init__(self, option_strings: Sequence[str], dest: str = ..., default: str = ..., help: str | None = ...) -> None: ...

# undocumented
class _VersionAction(Action):
    version: str | None
    def __init__(
        self, option_strings: Sequence[str], version: str | None = ..., dest: str = ..., default: str = ..., help: str = ...
    ) -> None: ...

# undocumented
class _SubParsersAction(Action, Generic[_ArgumentParserT]):
    _ChoicesPseudoAction: Type[Any]  # nested class
    _prog_prefix: str
    _parser_class: Type[_ArgumentParserT]
    _name_parser_map: dict[str, _ArgumentParserT]
    choices: dict[str, _ArgumentParserT]
    _choices_actions: list[Action]
    if sys.version_info >= (3, 7):
        def __init__(
            self,
            option_strings: Sequence[str],
            prog: str,
            parser_class: Type[_ArgumentParserT],
            dest: str = ...,
            required: bool = ...,
            help: str | None = ...,
            metavar: str | Tuple[str, ...] | None = ...,
        ) -> None: ...
    else:
        def __init__(
            self,
            option_strings: Sequence[str],
            prog: str,
            parser_class: Type[_ArgumentParserT],
            dest: str = ...,
            help: str | None = ...,
            metavar: str | Tuple[str, ...] | None = ...,
        ) -> None: ...
    # TODO: Type keyword args properly.
    def add_parser(self, name: str, **kwargs: Any) -> _ArgumentParserT: ...
    def _get_subactions(self) -> list[Action]: ...

# undocumented
class ArgumentTypeError(Exception): ...

if sys.version_info < (3, 7):
    # undocumented
    def _ensure_value(namespace: Namespace, name: str, value: Any) -> Any: ...

# undocumented
def _get_action_name(argument: Action | None) -> str | None: ...
