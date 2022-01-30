import os.path
import sys
import traceback
from mypy.backports import OrderedDict
from collections import defaultdict

from typing import Tuple, List, TypeVar, Set, Dict, Optional, TextIO, Callable
from typing_extensions import Final

from mypy.scope import Scope
from mypy.options import Options
from mypy.version import __version__ as mypy_version
from mypy.errorcodes import ErrorCode, IMPORT
from mypy.message_registry import ErrorMessage
from mypy import errorcodes as codes
from mypy.util import DEFAULT_SOURCE_OFFSET, is_typeshed_file

T = TypeVar("T")
allowed_duplicates: Final = ["@overload", "Got:", "Expected:"]


class ErrorInfo:
    """Representation of a single error message."""

    # Description of a sequence of imports that refer to the source file
    # related to this error. Each item is a (path, line number) tuple.
    import_ctx: List[Tuple[str, int]]

    # The path to source file that was the source of this error.
    file = ''

    # The fully-qualified id of the source module for this error.
    module: Optional[str] = None

    # The name of the type in which this error is located at.
    type: Optional[str] = ""  # Unqualified, may be None

    # The name of the function or member in which this error is located at.
    function_or_member: Optional[str] = ""  # Unqualified, may be None

    # The line number related to this error within file.
    line = 0     # -1 if unknown

    # The column number related to this error with file.
    column = 0   # -1 if unknown

    # Either 'error' or 'note'
    severity = ''

    # The error message.
    message = ''

    # The error code.
    code: Optional[ErrorCode] = None

    # If True, we should halt build after the file that generated this error.
    blocker = False

    # Only report this particular messages once per program.
    only_once = False

    # Do not remove duplicate copies of this message (ignored if only_once is True).
    allow_dups = False

    # Actual origin of the error message as tuple (path, line number, end line number)
    # If end line number is unknown, use line number.
    origin: Tuple[str, int, int]

    # Fine-grained incremental target where this was reported
    target: Optional[str] = None

    # If True, don't show this message in output, but still record the error (needed
    # by mypy daemon)
    hidden = False

    def __init__(self,
                 import_ctx: List[Tuple[str, int]],
                 file: str,
                 module: Optional[str],
                 typ: Optional[str],
                 function_or_member: Optional[str],
                 line: int,
                 column: int,
                 severity: str,
                 message: str,
                 code: Optional[ErrorCode],
                 blocker: bool,
                 only_once: bool,
                 allow_dups: bool,
                 origin: Optional[Tuple[str, int, int]] = None,
                 target: Optional[str] = None) -> None:
        self.import_ctx = import_ctx
        self.file = file
        self.module = module
        self.type = typ
        self.function_or_member = function_or_member
        self.line = line
        self.column = column
        self.severity = severity
        self.message = message
        self.code = code
        self.blocker = blocker
        self.only_once = only_once
        self.allow_dups = allow_dups
        self.origin = origin or (file, line, line)
        self.target = target


# Type used internally to represent errors:
#   (path, line, column, severity, message, allow_dups, code)
ErrorTuple = Tuple[Optional[str],
                   int,
                   int,
                   str,
                   str,
                   bool,
                   Optional[ErrorCode]]


class Errors:
    """Container for compile errors.

    This class generates and keeps tracks of compile errors and the
    current error context (nested imports).
    """

    # Map from files to generated error messages. Is an OrderedDict so
    # that it can be used to order messages based on the order the
    # files were processed.
    error_info_map: Dict[str, List[ErrorInfo]]

    # Files that we have reported the errors for
    flushed_files: Set[str]

    # Current error context: nested import context/stack, as a list of (path, line) pairs.
    import_ctx: List[Tuple[str, int]]

    # Path name prefix that is removed from all paths, if set.
    ignore_prefix: Optional[str] = None

    # Path to current file.
    file: str = ""

    # Ignore some errors on these lines of each file
    # (path -> line -> error-codes)
    ignored_lines: Dict[str, Dict[int, List[str]]]

    # Lines on which an error was actually ignored.
    used_ignored_lines: Dict[str, Dict[int, List[str]]]

    # Files where all errors should be ignored.
    ignored_files: Set[str]

    # Collection of reported only_once messages.
    only_once_messages: Set[str]

    # Set to True to show "In function "foo":" messages.
    show_error_context: bool = False

    # Set to True to show column numbers in error messages.
    show_column_numbers: bool = False

    # Set to True to show absolute file paths in error messages.
    show_absolute_path: bool = False

    # State for keeping track of the current fine-grained incremental mode target.
    # (See mypy.server.update for more about targets.)
    # Current module id.
    target_module: Optional[str] = None
    scope: Optional[Scope] = None

    # Have we seen an import-related error so far? If yes, we filter out other messages
    # in some cases to avoid reporting huge numbers of errors.
    seen_import_error = False

    def __init__(self,
                 show_error_context: bool = False,
                 show_column_numbers: bool = False,
                 show_error_codes: bool = False,
                 pretty: bool = False,
                 read_source: Optional[Callable[[str], Optional[List[str]]]] = None,
                 show_absolute_path: bool = False,
                 enabled_error_codes: Optional[Set[ErrorCode]] = None,
                 disabled_error_codes: Optional[Set[ErrorCode]] = None,
                 many_errors_threshold: int = -1) -> None:
        self.show_error_context = show_error_context
        self.show_column_numbers = show_column_numbers
        self.show_error_codes = show_error_codes
        self.show_absolute_path = show_absolute_path
        self.pretty = pretty
        # We use fscache to read source code when showing snippets.
        self.read_source = read_source
        self.enabled_error_codes = enabled_error_codes or set()
        self.disabled_error_codes = disabled_error_codes or set()
        self.many_errors_threshold = many_errors_threshold
        self.initialize()

    def initialize(self) -> None:
        self.error_info_map = OrderedDict()
        self.flushed_files = set()
        self.import_ctx = []
        self.function_or_member = [None]
        self.ignored_lines = OrderedDict()
        self.used_ignored_lines = defaultdict(lambda: defaultdict(list))
        self.ignored_files = set()
        self.only_once_messages = set()
        self.scope = None
        self.target_module = None
        self.seen_import_error = False

    def reset(self) -> None:
        self.initialize()

    def copy(self) -> 'Errors':
        new = Errors(self.show_error_context,
                     self.show_column_numbers,
                     self.show_error_codes,
                     self.pretty,
                     self.read_source,
                     self.show_absolute_path,
                     self.enabled_error_codes,
                     self.disabled_error_codes,
                     self.many_errors_threshold)
        new.file = self.file
        new.import_ctx = self.import_ctx[:]
        new.function_or_member = self.function_or_member[:]
        new.target_module = self.target_module
        new.scope = self.scope
        new.seen_import_error = self.seen_import_error
        return new

    def total_errors(self) -> int:
        return sum(len(errs) for errs in self.error_info_map.values())

    def set_ignore_prefix(self, prefix: str) -> None:
        """Set path prefix that will be removed from all paths."""
        prefix = os.path.normpath(prefix)
        # Add separator to the end, if not given.
        if os.path.basename(prefix) != '':
            prefix += os.sep
        self.ignore_prefix = prefix

    def simplify_path(self, file: str) -> str:
        if self.show_absolute_path:
            return os.path.abspath(file)
        else:
            file = os.path.normpath(file)
            return remove_path_prefix(file, self.ignore_prefix)

    def set_file(self, file: str,
                 module: Optional[str],
                 scope: Optional[Scope] = None) -> None:
        """Set the path and module id of the current file."""
        # The path will be simplified later, in render_messages. That way
        #  * 'file' is always a key that uniquely identifies a source file
        #    that mypy read (simplified paths might not be unique); and
        #  * we only have to simplify in one place, while still supporting
        #    reporting errors for files other than the one currently being
        #    processed.
        self.file = file
        self.target_module = module
        self.scope = scope

    def set_file_ignored_lines(self, file: str,
                               ignored_lines: Dict[int, List[str]],
                               ignore_all: bool = False) -> None:
        self.ignored_lines[file] = ignored_lines
        if ignore_all:
            self.ignored_files.add(file)

    def current_target(self) -> Optional[str]:
        """Retrieves the current target from the associated scope.

        If there is no associated scope, use the target module."""
        if self.scope is not None:
            return self.scope.current_target()
        return self.target_module

    def current_module(self) -> Optional[str]:
        return self.target_module

    def import_context(self) -> List[Tuple[str, int]]:
        """Return a copy of the import context."""
        return self.import_ctx[:]

    def set_import_context(self, ctx: List[Tuple[str, int]]) -> None:
        """Replace the entire import context with a new value."""
        self.import_ctx = ctx[:]

    def report(self,
               line: int,
               column: Optional[int],
               message: str,
               code: Optional[ErrorCode] = None,
               *,
               blocker: bool = False,
               severity: str = 'error',
               file: Optional[str] = None,
               only_once: bool = False,
               allow_dups: bool = False,
               origin_line: Optional[int] = None,
               offset: int = 0,
               end_line: Optional[int] = None) -> None:
        """Report message at the given line using the current error context.

        Args:
            line: line number of error
            column: column number of error
            message: message to report
            code: error code (defaults to 'misc'; not shown for notes)
            blocker: if True, don't continue analysis after this error
            severity: 'error' or 'note'
            file: if non-None, override current file as context
            only_once: if True, only report this exact message once per build
            allow_dups: if True, allow duplicate copies of this message (ignored if only_once)
            origin_line: if non-None, override current context as origin
            end_line: if non-None, override current context as end
        """
        if self.scope:
            type = self.scope.current_type_name()
            if self.scope.ignored > 0:
                type = None  # Omit type context if nested function
            function = self.scope.current_function_name()
        else:
            type = None
            function = None

        if column is None:
            column = -1
        if file is None:
            file = self.file
        if offset:
            message = " " * offset + message

        if origin_line is None:
            origin_line = line

        if end_line is None:
            end_line = origin_line

        code = code or (codes.MISC if not blocker else None)

        info = ErrorInfo(self.import_context(), file, self.current_module(), type,
                         function, line, column, severity, message, code,
                         blocker, only_once, allow_dups,
                         origin=(self.file, origin_line, end_line),
                         target=self.current_target())
        self.add_error_info(info)

    def _add_error_info(self, file: str, info: ErrorInfo) -> None:
        assert file not in self.flushed_files
        if file not in self.error_info_map:
            self.error_info_map[file] = []
        self.error_info_map[file].append(info)
        if info.code is IMPORT:
            self.seen_import_error = True

    def add_error_info(self, info: ErrorInfo) -> None:
        file, line, end_line = info.origin
        if not info.blocker:  # Blockers cannot be ignored
            if file in self.ignored_lines:
                # It's okay if end_line is *before* line.
                # Function definitions do this, for example, because the correct
                # error reporting line is at the *end* of the ignorable range
                # (for compatibility reasons). If so, just flip 'em!
                if end_line < line:
                    line, end_line = end_line, line
                # Check each line in this context for "type: ignore" comments.
                # line == end_line for most nodes, so we only loop once.
                for scope_line in range(line, end_line + 1):
                    if self.is_ignored_error(scope_line, info, self.ignored_lines[file]):
                        # Annotation requests us to ignore all errors on this line.
                        self.used_ignored_lines[file][scope_line].append(
                            (info.code or codes.MISC).code)
                        return
            if file in self.ignored_files:
                return
        if info.only_once:
            if info.message in self.only_once_messages:
                return
            self.only_once_messages.add(info.message)
        if self.seen_import_error and info.code is not IMPORT and self.has_many_errors():
            # Missing stubs can easily cause thousands of errors about
            # Any types, especially when upgrading to mypy 0.900,
            # which no longer bundles third-party library stubs. Avoid
            # showing too many errors to make it easier to see
            # import-related errors.
            info.hidden = True
            self.report_hidden_errors(info)
        self._add_error_info(file, info)

    def has_many_errors(self) -> bool:
        if self.many_errors_threshold < 0:
            return False
        if len(self.error_info_map) >= self.many_errors_threshold:
            return True
        if sum(len(errors)
               for errors in self.error_info_map.values()) >= self.many_errors_threshold:
            return True
        return False

    def report_hidden_errors(self, info: ErrorInfo) -> None:
        message = (
            '(Skipping most remaining errors due to unresolved imports or missing stubs; ' +
            'fix these first)'
        )
        if message in self.only_once_messages:
            return
        self.only_once_messages.add(message)
        new_info = ErrorInfo(
            import_ctx=info.import_ctx,
            file=info.file,
            module=info.module,
            typ=None,
            function_or_member=None,
            line=info.line,
            column=info.line,
            severity='note',
            message=message,
            code=None,
            blocker=False,
            only_once=True,
            allow_dups=False,
            origin=info.origin,
            target=info.target,
        )
        self._add_error_info(info.origin[0], new_info)

    def is_ignored_error(self, line: int, info: ErrorInfo, ignores: Dict[int, List[str]]) -> bool:
        if info.blocker:
            # Blocking errors can never be ignored
            return False
        if info.code and self.is_error_code_enabled(info.code) is False:
            return True
        if line not in ignores:
            return False
        if not ignores[line]:
            # Empty list means that we ignore all errors
            return True
        if info.code and self.is_error_code_enabled(info.code) is True:
            return info.code.code in ignores[line]
        return False

    def is_error_code_enabled(self, error_code: ErrorCode) -> bool:
        if error_code in self.disabled_error_codes:
            return False
        elif error_code in self.enabled_error_codes:
            return True
        else:
            return error_code.default_enabled

    def clear_errors_in_targets(self, path: str, targets: Set[str]) -> None:
        """Remove errors in specific fine-grained targets within a file."""
        if path in self.error_info_map:
            new_errors = []
            for info in self.error_info_map[path]:
                if info.target not in targets:
                    new_errors.append(info)
                elif info.only_once:
                    self.only_once_messages.remove(info.message)
            self.error_info_map[path] = new_errors

    def generate_unused_ignore_errors(self, file: str) -> None:
        ignored_lines = self.ignored_lines[file]
        if not is_typeshed_file(file) and file not in self.ignored_files:
            ignored_lines = self.ignored_lines[file]
            used_ignored_lines = self.used_ignored_lines[file]
            for line, ignored_codes in ignored_lines.items():
                used_ignored_codes = used_ignored_lines[line]
                unused_ignored_codes = set(ignored_codes) - set(used_ignored_codes)
                # `ignore` is used
                if len(ignored_codes) == 0 and len(used_ignored_codes) > 0:
                    continue
                # All codes appearing in `ignore[...]` are used
                if len(ignored_codes) > 0 and len(unused_ignored_codes) == 0:
                    continue
                # Display detail only when `ignore[...]` specifies more than one error code
                unused_codes_message = ""
                if len(ignored_codes) > 1 and len(unused_ignored_codes) > 0:
                    unused_codes_message = f"[{', '.join(sorted(unused_ignored_codes))}]"
                message = f'Unused "type: ignore{unused_codes_message}" comment'
                # Don't use report since add_error_info will ignore the error!
                info = ErrorInfo(self.import_context(), file, self.current_module(), None,
                                 None, line, -1, 'error', message,
                                 None, False, False, False)
                self._add_error_info(file, info)

    def num_messages(self) -> int:
        """Return the number of generated messages."""
        return sum(len(x) for x in self.error_info_map.values())

    def is_errors(self) -> bool:
        """Are there any generated messages?"""
        return bool(self.error_info_map)

    def is_real_errors(self) -> bool:
        """Are there any generated errors (not just notes, for example)?"""
        return any(info.severity == 'error'
                   for infos in self.error_info_map.values() for info in infos)

    def is_blockers(self) -> bool:
        """Are the any errors that are blockers?"""
        return any(err for errs in self.error_info_map.values() for err in errs if err.blocker)

    def blocker_module(self) -> Optional[str]:
        """Return the module with a blocking error, or None if not possible."""
        for errs in self.error_info_map.values():
            for err in errs:
                if err.blocker:
                    return err.module
        return None

    def is_errors_for_file(self, file: str) -> bool:
        """Are there any errors for the given file?"""
        return file in self.error_info_map

    def most_recent_error_location(self) -> Tuple[int, int]:
        info = self.error_info_map[self.file][-1]
        return info.line, info.column

    def raise_error(self, use_stdout: bool = True) -> None:
        """Raise a CompileError with the generated messages.

        Render the messages suitable for displaying.
        """
        # self.new_messages() will format all messages that haven't already
        # been returned from a file_messages() call.
        raise CompileError(self.new_messages(),
                           use_stdout=use_stdout,
                           module_with_blocker=self.blocker_module())

    def format_messages(self, error_info: List[ErrorInfo],
                        source_lines: Optional[List[str]]) -> List[str]:
        """Return a string list that represents the error messages.

        Use a form suitable for displaying to the user. If self.pretty
        is True also append a relevant trimmed source code line (only for
        severity 'error').
        """
        a: List[str] = []
        error_info = [info for info in error_info if not info.hidden]
        errors = self.render_messages(self.sort_messages(error_info))
        errors = self.remove_duplicates(errors)
        for file, line, column, severity, message, allow_dups, code in errors:
            s = ''
            if file is not None:
                if self.show_column_numbers and line >= 0 and column >= 0:
                    srcloc = '{}:{}:{}'.format(file, line, 1 + column)
                elif line >= 0:
                    srcloc = '{}:{}'.format(file, line)
                else:
                    srcloc = file
                s = '{}: {}: {}'.format(srcloc, severity, message)
            else:
                s = message
            if self.show_error_codes and code and severity != 'note':
                # If note has an error code, it is related to a previous error. Avoid
                # displaying duplicate error codes.
                s = '{}  [{}]'.format(s, code.code)
            a.append(s)
            if self.pretty:
                # Add source code fragment and a location marker.
                if severity == 'error' and source_lines and line > 0:
                    source_line = source_lines[line - 1]
                    source_line_expanded = source_line.expandtabs()
                    if column < 0:
                        # Something went wrong, take first non-empty column.
                        column = len(source_line) - len(source_line.lstrip())

                    # Shifts column after tab expansion
                    column = len(source_line[:column].expandtabs())

                    # Note, currently coloring uses the offset to detect source snippets,
                    # so these offsets should not be arbitrary.
                    a.append(' ' * DEFAULT_SOURCE_OFFSET + source_line_expanded)
                    a.append(' ' * (DEFAULT_SOURCE_OFFSET + column) + '^')
        return a

    def file_messages(self, path: str) -> List[str]:
        """Return a string list of new error messages from a given file.

        Use a form suitable for displaying to the user.
        """
        if path not in self.error_info_map:
            return []
        self.flushed_files.add(path)
        source_lines = None
        if self.pretty:
            assert self.read_source
            source_lines = self.read_source(path)
        return self.format_messages(self.error_info_map[path], source_lines)

    def new_messages(self) -> List[str]:
        """Return a string list of new error messages.

        Use a form suitable for displaying to the user.
        Errors from different files are ordered based on the order in which
        they first generated an error.
        """
        msgs = []
        for path in self.error_info_map.keys():
            if path not in self.flushed_files:
                msgs.extend(self.file_messages(path))
        return msgs

    def targets(self) -> Set[str]:
        """Return a set of all targets that contain errors."""
        # TODO: Make sure that either target is always defined or that not being defined
        #       is okay for fine-grained incremental checking.
        return set(info.target
                   for errs in self.error_info_map.values()
                   for info in errs
                   if info.target)

    def render_messages(self,
                        errors: List[ErrorInfo]) -> List[ErrorTuple]:
        """Translate the messages into a sequence of tuples.

        Each tuple is of form (path, line, col, severity, message, allow_dups, code).
        The rendered sequence includes information about error contexts.
        The path item may be None. If the line item is negative, the
        line number is not defined for the tuple.
        """
        result: List[ErrorTuple] = []
        prev_import_context: List[Tuple[str, int]] = []
        prev_function_or_member: Optional[str] = None
        prev_type: Optional[str] = None

        for e in errors:
            # Report module import context, if different from previous message.
            if not self.show_error_context:
                pass
            elif e.import_ctx != prev_import_context:
                last = len(e.import_ctx) - 1
                i = last
                while i >= 0:
                    path, line = e.import_ctx[i]
                    fmt = '{}:{}: note: In module imported here'
                    if i < last:
                        fmt = '{}:{}: note: ... from here'
                    if i > 0:
                        fmt += ','
                    else:
                        fmt += ':'
                    # Remove prefix to ignore from path (if present) to
                    # simplify path.
                    path = remove_path_prefix(path, self.ignore_prefix)
                    result.append((None, -1, -1, 'note',
                                   fmt.format(path, line), e.allow_dups, None))
                    i -= 1

            file = self.simplify_path(e.file)

            # Report context within a source file.
            if not self.show_error_context:
                pass
            elif (e.function_or_member != prev_function_or_member or
                    e.type != prev_type):
                if e.function_or_member is None:
                    if e.type is None:
                        result.append((file, -1, -1, 'note', 'At top level:', e.allow_dups, None))
                    else:
                        result.append((file, -1, -1, 'note', 'In class "{}":'.format(
                            e.type), e.allow_dups, None))
                else:
                    if e.type is None:
                        result.append((file, -1, -1, 'note',
                                       'In function "{}":'.format(
                                           e.function_or_member), e.allow_dups, None))
                    else:
                        result.append((file, -1, -1, 'note',
                                       'In member "{}" of class "{}":'.format(
                                           e.function_or_member, e.type), e.allow_dups, None))
            elif e.type != prev_type:
                if e.type is None:
                    result.append((file, -1, -1, 'note', 'At top level:', e.allow_dups, None))
                else:
                    result.append((file, -1, -1, 'note',
                                   'In class "{}":'.format(e.type), e.allow_dups, None))

            if isinstance(e.message, ErrorMessage):
                result.append(
                    (file, e.line, e.column, e.severity, e.message.value, e.allow_dups, e.code))
            else:
                result.append(
                    (file, e.line, e.column, e.severity, e.message, e.allow_dups, e.code))

            prev_import_context = e.import_ctx
            prev_function_or_member = e.function_or_member
            prev_type = e.type

        return result

    def sort_messages(self, errors: List[ErrorInfo]) -> List[ErrorInfo]:
        """Sort an array of error messages locally by line number.

        I.e., sort a run of consecutive messages with the same
        context by line number, but otherwise retain the general
        ordering of the messages.
        """
        result: List[ErrorInfo] = []
        i = 0
        while i < len(errors):
            i0 = i
            # Find neighbouring errors with the same context and file.
            while (i + 1 < len(errors) and
                    errors[i + 1].import_ctx == errors[i].import_ctx and
                    errors[i + 1].file == errors[i].file):
                i += 1
            i += 1

            # Sort the errors specific to a file according to line number and column.
            a = sorted(errors[i0:i], key=lambda x: (x.line, x.column))
            result.extend(a)
        return result

    def remove_duplicates(self, errors: List[ErrorTuple]) -> List[ErrorTuple]:
        """Remove duplicates from a sorted error list."""
        res: List[ErrorTuple] = []
        i = 0
        while i < len(errors):
            dup = False
            # Use slightly special formatting for member conflicts reporting.
            conflicts_notes = False
            j = i - 1
            # Find duplicates, unless duplicates are allowed.
            if not errors[i][5]:
                while j >= 0 and errors[j][0] == errors[i][0]:
                    if errors[j][4].strip() == 'Got:':
                        conflicts_notes = True
                    j -= 1
                j = i - 1
                while (j >= 0 and errors[j][0] == errors[i][0] and
                        errors[j][1] == errors[i][1]):
                    if (errors[j][3] == errors[i][3] and
                            # Allow duplicate notes in overload conflicts reporting.
                            not ((errors[i][3] == 'note' and
                                errors[i][4].strip() in allowed_duplicates)
                                or (errors[i][4].strip().startswith('def ') and
                                    conflicts_notes)) and
                            errors[j][4] == errors[i][4]):  # ignore column
                        dup = True
                        break
                    j -= 1
            if not dup:
                res.append(errors[i])
            i += 1
        return res


class CompileError(Exception):
    """Exception raised when there is a compile error.

    It can be a parse, semantic analysis, type check or other
    compilation-related error.

    CompileErrors raised from an errors object carry all of the
    messages that have not been reported out by error streaming.
    This is patched up by build.build to contain either all error
    messages (if errors were streamed) or none (if they were not).

    """

    messages: List[str]
    use_stdout = False
    # Can be set in case there was a module with a blocking error
    module_with_blocker: Optional[str] = None

    def __init__(self,
                 messages: List[str],
                 use_stdout: bool = False,
                 module_with_blocker: Optional[str] = None) -> None:
        super().__init__('\n'.join(messages))
        self.messages = messages
        self.use_stdout = use_stdout
        self.module_with_blocker = module_with_blocker


def remove_path_prefix(path: str, prefix: Optional[str]) -> str:
    """If path starts with prefix, return copy of path with the prefix removed.
    Otherwise, return path. If path is None, return None.
    """
    if prefix is not None and path.startswith(prefix):
        return path[len(prefix):]
    else:
        return path


def report_internal_error(err: Exception,
                          file: Optional[str],
                          line: int,
                          errors: Errors,
                          options: Options,
                          stdout: Optional[TextIO] = None,
                          stderr: Optional[TextIO] = None,
                          ) -> None:
    """Report internal error and exit.

    This optionally starts pdb or shows a traceback.
    """
    stdout = (stdout or sys.stdout)
    stderr = (stderr or sys.stderr)
    # Dump out errors so far, they often provide a clue.
    # But catch unexpected errors rendering them.
    try:
        for msg in errors.new_messages():
            print(msg)
    except Exception as e:
        print("Failed to dump errors:", repr(e), file=stderr)

    # Compute file:line prefix for official-looking error messages.
    if file:
        if line:
            prefix = '{}:{}: '.format(file, line)
        else:
            prefix = '{}: '.format(file)
    else:
        prefix = ''

    # Print "INTERNAL ERROR" message.
    print('{}error: INTERNAL ERROR --'.format(prefix),
          'Please try using mypy master on Github:\n'
          'https://mypy.readthedocs.io/en/stable/common_issues.html'
          '#using-a-development-mypy-build',
          file=stderr)
    if options.show_traceback:
        print('Please report a bug at https://github.com/python/mypy/issues',
            file=stderr)
    else:
        print('If this issue continues with mypy master, '
              'please report a bug at https://github.com/python/mypy/issues',
            file=stderr)
    print('version: {}'.format(mypy_version),
          file=stderr)

    # If requested, drop into pdb. This overrides show_tb.
    if options.pdb:
        print('Dropping into pdb', file=stderr)
        import pdb
        pdb.post_mortem(sys.exc_info()[2])

    # If requested, print traceback, else print note explaining how to get one.
    if options.raise_exceptions:
        raise err
    if not options.show_traceback:
        if not options.pdb:
            print('{}: note: please use --show-traceback to print a traceback '
                  'when reporting a bug'.format(prefix),
                  file=stderr)
    else:
        tb = traceback.extract_stack()[:-2]
        tb2 = traceback.extract_tb(sys.exc_info()[2])
        print('Traceback (most recent call last):')
        for s in traceback.format_list(tb + tb2):
            print(s.rstrip('\n'))
        print('{}: {}'.format(type(err).__name__, err), file=stdout)
        print('{}: note: use --pdb to drop into pdb'.format(prefix), file=stderr)

    # Exit.  The caller has nothing more to say.
    # We use exit code 2 to signal that this is no ordinary error.
    raise SystemExit(2)
