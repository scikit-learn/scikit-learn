import asyncio
from concurrent.futures import Executor, ThreadPoolExecutor, ProcessPoolExecutor
from contextlib import contextmanager
from datetime import datetime
from enum import Enum
import io
from multiprocessing import Manager, freeze_support
import os
from pathlib import Path
import regex as re
import signal
import sys
import tokenize
import traceback
from typing import (
    Any,
    Dict,
    Generator,
    Iterator,
    List,
    Optional,
    Pattern,
    Set,
    Sized,
    Tuple,
    Union,
)

from dataclasses import replace
import click

from black.const import DEFAULT_LINE_LENGTH, DEFAULT_INCLUDES, DEFAULT_EXCLUDES
from black.const import STDIN_PLACEHOLDER
from black.nodes import STARS, syms, is_simple_decorator_expression
from black.lines import Line, EmptyLineTracker
from black.linegen import transform_line, LineGenerator, LN
from black.comments import normalize_fmt_off
from black.mode import Mode, TargetVersion
from black.mode import Feature, supports_feature, VERSION_TO_FEATURES
from black.cache import read_cache, write_cache, get_cache_info, filter_cached, Cache
from black.concurrency import cancel, shutdown, maybe_install_uvloop
from black.output import dump_to_file, diff, color_diff, out, err
from black.report import Report, Changed
from black.files import find_project_root, find_pyproject_toml, parse_pyproject_toml
from black.files import gen_python_files, get_gitignore, normalize_path_maybe_ignore
from black.files import wrap_stream_for_windows
from black.parsing import InvalidInput  # noqa F401
from black.parsing import lib2to3_parse, parse_ast, stringify_ast


# lib2to3 fork
from blib2to3.pytree import Node, Leaf
from blib2to3.pgen2 import token

from _black_version import version as __version__

# types
FileContent = str
Encoding = str
NewLine = str


class NothingChanged(UserWarning):
    """Raised when reformatted code is the same as source."""


class WriteBack(Enum):
    NO = 0
    YES = 1
    DIFF = 2
    CHECK = 3
    COLOR_DIFF = 4

    @classmethod
    def from_configuration(
        cls, *, check: bool, diff: bool, color: bool = False
    ) -> "WriteBack":
        if check and not diff:
            return cls.CHECK

        if diff and color:
            return cls.COLOR_DIFF

        return cls.DIFF if diff else cls.YES


# Legacy name, left for integrations.
FileMode = Mode


def read_pyproject_toml(
    ctx: click.Context, param: click.Parameter, value: Optional[str]
) -> Optional[str]:
    """Inject Black configuration from "pyproject.toml" into defaults in `ctx`.

    Returns the path to a successfully found and read configuration file, None
    otherwise.
    """
    if not value:
        value = find_pyproject_toml(ctx.params.get("src", ()))
        if value is None:
            return None

    try:
        config = parse_pyproject_toml(value)
    except (OSError, ValueError) as e:
        raise click.FileError(
            filename=value, hint=f"Error reading configuration file: {e}"
        )

    if not config:
        return None
    else:
        # Sanitize the values to be Click friendly. For more information please see:
        # https://github.com/psf/black/issues/1458
        # https://github.com/pallets/click/issues/1567
        config = {
            k: str(v) if not isinstance(v, (list, dict)) else v
            for k, v in config.items()
        }

    target_version = config.get("target_version")
    if target_version is not None and not isinstance(target_version, list):
        raise click.BadOptionUsage(
            "target-version", "Config key target-version must be a list"
        )

    default_map: Dict[str, Any] = {}
    if ctx.default_map:
        default_map.update(ctx.default_map)
    default_map.update(config)

    ctx.default_map = default_map
    return value


def target_version_option_callback(
    c: click.Context, p: Union[click.Option, click.Parameter], v: Tuple[str, ...]
) -> List[TargetVersion]:
    """Compute the target versions from a --target-version flag.

    This is its own function because mypy couldn't infer the type correctly
    when it was a lambda, causing mypyc trouble.
    """
    return [TargetVersion[val.upper()] for val in v]


def re_compile_maybe_verbose(regex: str) -> Pattern[str]:
    """Compile a regular expression string in `regex`.

    If it contains newlines, use verbose mode.
    """
    if "\n" in regex:
        regex = "(?x)" + regex
    compiled: Pattern[str] = re.compile(regex)
    return compiled


def validate_regex(
    ctx: click.Context,
    param: click.Parameter,
    value: Optional[str],
) -> Optional[Pattern]:
    try:
        return re_compile_maybe_verbose(value) if value is not None else None
    except re.error:
        raise click.BadParameter("Not a valid regular expression")


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option("-c", "--code", type=str, help="Format the code passed in as a string.")
@click.option(
    "-l",
    "--line-length",
    type=int,
    default=DEFAULT_LINE_LENGTH,
    help="How many characters per line to allow.",
    show_default=True,
)
@click.option(
    "-t",
    "--target-version",
    type=click.Choice([v.name.lower() for v in TargetVersion]),
    callback=target_version_option_callback,
    multiple=True,
    help=(
        "Python versions that should be supported by Black's output. [default: per-file"
        " auto-detection]"
    ),
)
@click.option(
    "--pyi",
    is_flag=True,
    help=(
        "Format all input files like typing stubs regardless of file extension (useful"
        " when piping source on standard input)."
    ),
)
@click.option(
    "-S",
    "--skip-string-normalization",
    is_flag=True,
    help="Don't normalize string quotes or prefixes.",
)
@click.option(
    "-C",
    "--skip-magic-trailing-comma",
    is_flag=True,
    help="Don't use trailing commas as a reason to split lines.",
)
@click.option(
    "--experimental-string-processing",
    is_flag=True,
    hidden=True,
    help=(
        "Experimental option that performs more normalization on string literals."
        " Currently disabled because it leads to some crashes."
    ),
)
@click.option(
    "--check",
    is_flag=True,
    help=(
        "Don't write the files back, just return the status. Return code 0 means"
        " nothing would change. Return code 1 means some files would be reformatted."
        " Return code 123 means there was an internal error."
    ),
)
@click.option(
    "--diff",
    is_flag=True,
    help="Don't write the files back, just output a diff for each file on stdout.",
)
@click.option(
    "--color/--no-color",
    is_flag=True,
    help="Show colored diff. Only applies when `--diff` is given.",
)
@click.option(
    "--fast/--safe",
    is_flag=True,
    help="If --fast given, skip temporary sanity checks. [default: --safe]",
)
@click.option(
    "--required-version",
    type=str,
    help=(
        "Require a specific version of Black to be running (useful for unifying results"
        " across many environments e.g. with a pyproject.toml file)."
    ),
)
@click.option(
    "--include",
    type=str,
    default=DEFAULT_INCLUDES,
    callback=validate_regex,
    help=(
        "A regular expression that matches files and directories that should be"
        " included on recursive searches. An empty value means all files are included"
        " regardless of the name. Use forward slashes for directories on all platforms"
        " (Windows, too). Exclusions are calculated first, inclusions later."
    ),
    show_default=True,
)
@click.option(
    "--exclude",
    type=str,
    callback=validate_regex,
    help=(
        "A regular expression that matches files and directories that should be"
        " excluded on recursive searches. An empty value means no paths are excluded."
        " Use forward slashes for directories on all platforms (Windows, too)."
        " Exclusions are calculated first, inclusions later. [default:"
        f" {DEFAULT_EXCLUDES}]"
    ),
    show_default=False,
)
@click.option(
    "--extend-exclude",
    type=str,
    callback=validate_regex,
    help=(
        "Like --exclude, but adds additional files and directories on top of the"
        " excluded ones. (Useful if you simply want to add to the default)"
    ),
)
@click.option(
    "--force-exclude",
    type=str,
    callback=validate_regex,
    help=(
        "Like --exclude, but files and directories matching this regex will be "
        "excluded even when they are passed explicitly as arguments."
    ),
)
@click.option(
    "--stdin-filename",
    type=str,
    help=(
        "The name of the file when passing it through stdin. Useful to make "
        "sure Black will respect --force-exclude option on some "
        "editors that rely on using stdin."
    ),
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help=(
        "Don't emit non-error messages to stderr. Errors are still emitted; silence"
        " those with 2>/dev/null."
    ),
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help=(
        "Also emit messages to stderr about files that were not changed or were ignored"
        " due to exclusion patterns."
    ),
)
@click.version_option(version=__version__)
@click.argument(
    "src",
    nargs=-1,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=True, readable=True, allow_dash=True
    ),
    is_eager=True,
)
@click.option(
    "--config",
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        allow_dash=False,
        path_type=str,
    ),
    is_eager=True,
    callback=read_pyproject_toml,
    help="Read configuration from FILE path.",
)
@click.pass_context
def main(
    ctx: click.Context,
    code: Optional[str],
    line_length: int,
    target_version: List[TargetVersion],
    check: bool,
    diff: bool,
    color: bool,
    fast: bool,
    pyi: bool,
    skip_string_normalization: bool,
    skip_magic_trailing_comma: bool,
    experimental_string_processing: bool,
    quiet: bool,
    verbose: bool,
    required_version: str,
    include: Pattern,
    exclude: Optional[Pattern],
    extend_exclude: Optional[Pattern],
    force_exclude: Optional[Pattern],
    stdin_filename: Optional[str],
    src: Tuple[str, ...],
    config: Optional[str],
) -> None:
    """The uncompromising code formatter."""
    if config and verbose:
        out(f"Using configuration from {config}.", bold=False, fg="blue")

    error_msg = "Oh no! 💥 💔 💥"
    if required_version and required_version != __version__:
        err(
            f"{error_msg} The required version `{required_version}` does not match"
            f" the running version `{__version__}`!"
        )
        ctx.exit(1)

    write_back = WriteBack.from_configuration(check=check, diff=diff, color=color)
    if target_version:
        versions = set(target_version)
    else:
        # We'll autodetect later.
        versions = set()
    mode = Mode(
        target_versions=versions,
        line_length=line_length,
        is_pyi=pyi,
        string_normalization=not skip_string_normalization,
        magic_trailing_comma=not skip_magic_trailing_comma,
        experimental_string_processing=experimental_string_processing,
    )

    if code is not None:
        # Run in quiet mode by default with -c; the extra output isn't useful.
        # You can still pass -v to get verbose output.
        quiet = True

    report = Report(check=check, diff=diff, quiet=quiet, verbose=verbose)

    if code is not None:
        reformat_code(
            content=code, fast=fast, write_back=write_back, mode=mode, report=report
        )
    else:
        sources = get_sources(
            ctx=ctx,
            src=src,
            quiet=quiet,
            verbose=verbose,
            include=include,
            exclude=exclude,
            extend_exclude=extend_exclude,
            force_exclude=force_exclude,
            report=report,
            stdin_filename=stdin_filename,
        )

        path_empty(
            sources,
            "No Python files are present to be formatted. Nothing to do 😴",
            quiet,
            verbose,
            ctx,
        )

        if len(sources) == 1:
            reformat_one(
                src=sources.pop(),
                fast=fast,
                write_back=write_back,
                mode=mode,
                report=report,
            )
        else:
            reformat_many(
                sources=sources,
                fast=fast,
                write_back=write_back,
                mode=mode,
                report=report,
            )

    if verbose or not quiet:
        out(error_msg if report.return_code else "All done! ✨ 🍰 ✨")
        if code is None:
            click.echo(str(report), err=True)
    ctx.exit(report.return_code)


def get_sources(
    *,
    ctx: click.Context,
    src: Tuple[str, ...],
    quiet: bool,
    verbose: bool,
    include: Pattern[str],
    exclude: Optional[Pattern[str]],
    extend_exclude: Optional[Pattern[str]],
    force_exclude: Optional[Pattern[str]],
    report: "Report",
    stdin_filename: Optional[str],
) -> Set[Path]:
    """Compute the set of files to be formatted."""

    root = find_project_root(src)
    sources: Set[Path] = set()
    path_empty(src, "No Path provided. Nothing to do 😴", quiet, verbose, ctx)

    if exclude is None:
        exclude = re_compile_maybe_verbose(DEFAULT_EXCLUDES)
        gitignore = get_gitignore(root)
    else:
        gitignore = None

    for s in src:
        if s == "-" and stdin_filename:
            p = Path(stdin_filename)
            is_stdin = True
        else:
            p = Path(s)
            is_stdin = False

        if is_stdin or p.is_file():
            normalized_path = normalize_path_maybe_ignore(p, root, report)
            if normalized_path is None:
                continue

            normalized_path = "/" + normalized_path
            # Hard-exclude any files that matches the `--force-exclude` regex.
            if force_exclude:
                force_exclude_match = force_exclude.search(normalized_path)
            else:
                force_exclude_match = None
            if force_exclude_match and force_exclude_match.group(0):
                report.path_ignored(p, "matches the --force-exclude regular expression")
                continue

            if is_stdin:
                p = Path(f"{STDIN_PLACEHOLDER}{str(p)}")

            sources.add(p)
        elif p.is_dir():
            sources.update(
                gen_python_files(
                    p.iterdir(),
                    root,
                    include,
                    exclude,
                    extend_exclude,
                    force_exclude,
                    report,
                    gitignore,
                )
            )
        elif s == "-":
            sources.add(p)
        else:
            err(f"invalid path: {s}")
    return sources


def path_empty(
    src: Sized, msg: str, quiet: bool, verbose: bool, ctx: click.Context
) -> None:
    """
    Exit if there is no `src` provided for formatting
    """
    if not src:
        if verbose or not quiet:
            out(msg)
        ctx.exit(0)


def reformat_code(
    content: str, fast: bool, write_back: WriteBack, mode: Mode, report: Report
) -> None:
    """
    Reformat and print out `content` without spawning child processes.
    Similar to `reformat_one`, but for string content.

    `fast`, `write_back`, and `mode` options are passed to
    :func:`format_file_in_place` or :func:`format_stdin_to_stdout`.
    """
    path = Path("<string>")
    try:
        changed = Changed.NO
        if format_stdin_to_stdout(
            content=content, fast=fast, write_back=write_back, mode=mode
        ):
            changed = Changed.YES
        report.done(path, changed)
    except Exception as exc:
        if report.verbose:
            traceback.print_exc()
        report.failed(path, str(exc))


def reformat_one(
    src: Path, fast: bool, write_back: WriteBack, mode: Mode, report: "Report"
) -> None:
    """Reformat a single file under `src` without spawning child processes.

    `fast`, `write_back`, and `mode` options are passed to
    :func:`format_file_in_place` or :func:`format_stdin_to_stdout`.
    """
    try:
        changed = Changed.NO

        if str(src) == "-":
            is_stdin = True
        elif str(src).startswith(STDIN_PLACEHOLDER):
            is_stdin = True
            # Use the original name again in case we want to print something
            # to the user
            src = Path(str(src)[len(STDIN_PLACEHOLDER) :])
        else:
            is_stdin = False

        if is_stdin:
            if src.suffix == ".pyi":
                mode = replace(mode, is_pyi=True)
            if format_stdin_to_stdout(fast=fast, write_back=write_back, mode=mode):
                changed = Changed.YES
        else:
            cache: Cache = {}
            if write_back not in (WriteBack.DIFF, WriteBack.COLOR_DIFF):
                cache = read_cache(mode)
                res_src = src.resolve()
                res_src_s = str(res_src)
                if res_src_s in cache and cache[res_src_s] == get_cache_info(res_src):
                    changed = Changed.CACHED
            if changed is not Changed.CACHED and format_file_in_place(
                src, fast=fast, write_back=write_back, mode=mode
            ):
                changed = Changed.YES
            if (write_back is WriteBack.YES and changed is not Changed.CACHED) or (
                write_back is WriteBack.CHECK and changed is Changed.NO
            ):
                write_cache(cache, [src], mode)
        report.done(src, changed)
    except Exception as exc:
        if report.verbose:
            traceback.print_exc()
        report.failed(src, str(exc))


def reformat_many(
    sources: Set[Path], fast: bool, write_back: WriteBack, mode: Mode, report: "Report"
) -> None:
    """Reformat multiple files using a ProcessPoolExecutor."""
    executor: Executor
    loop = asyncio.get_event_loop()
    worker_count = os.cpu_count()
    if sys.platform == "win32":
        # Work around https://bugs.python.org/issue26903
        worker_count = min(worker_count, 60)
    try:
        executor = ProcessPoolExecutor(max_workers=worker_count)
    except (ImportError, OSError):
        # we arrive here if the underlying system does not support multi-processing
        # like in AWS Lambda or Termux, in which case we gracefully fallback to
        # a ThreadPoolExecutor with just a single worker (more workers would not do us
        # any good due to the Global Interpreter Lock)
        executor = ThreadPoolExecutor(max_workers=1)

    try:
        loop.run_until_complete(
            schedule_formatting(
                sources=sources,
                fast=fast,
                write_back=write_back,
                mode=mode,
                report=report,
                loop=loop,
                executor=executor,
            )
        )
    finally:
        shutdown(loop)
        if executor is not None:
            executor.shutdown()


async def schedule_formatting(
    sources: Set[Path],
    fast: bool,
    write_back: WriteBack,
    mode: Mode,
    report: "Report",
    loop: asyncio.AbstractEventLoop,
    executor: Executor,
) -> None:
    """Run formatting of `sources` in parallel using the provided `executor`.

    (Use ProcessPoolExecutors for actual parallelism.)

    `write_back`, `fast`, and `mode` options are passed to
    :func:`format_file_in_place`.
    """
    cache: Cache = {}
    if write_back not in (WriteBack.DIFF, WriteBack.COLOR_DIFF):
        cache = read_cache(mode)
        sources, cached = filter_cached(cache, sources)
        for src in sorted(cached):
            report.done(src, Changed.CACHED)
    if not sources:
        return

    cancelled = []
    sources_to_cache = []
    lock = None
    if write_back in (WriteBack.DIFF, WriteBack.COLOR_DIFF):
        # For diff output, we need locks to ensure we don't interleave output
        # from different processes.
        manager = Manager()
        lock = manager.Lock()
    tasks = {
        asyncio.ensure_future(
            loop.run_in_executor(
                executor, format_file_in_place, src, fast, mode, write_back, lock
            )
        ): src
        for src in sorted(sources)
    }
    pending = tasks.keys()
    try:
        loop.add_signal_handler(signal.SIGINT, cancel, pending)
        loop.add_signal_handler(signal.SIGTERM, cancel, pending)
    except NotImplementedError:
        # There are no good alternatives for these on Windows.
        pass
    while pending:
        done, _ = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
        for task in done:
            src = tasks.pop(task)
            if task.cancelled():
                cancelled.append(task)
            elif task.exception():
                report.failed(src, str(task.exception()))
            else:
                changed = Changed.YES if task.result() else Changed.NO
                # If the file was written back or was successfully checked as
                # well-formatted, store this information in the cache.
                if write_back is WriteBack.YES or (
                    write_back is WriteBack.CHECK and changed is Changed.NO
                ):
                    sources_to_cache.append(src)
                report.done(src, changed)
    if cancelled:
        await asyncio.gather(*cancelled, loop=loop, return_exceptions=True)
    if sources_to_cache:
        write_cache(cache, sources_to_cache, mode)


def format_file_in_place(
    src: Path,
    fast: bool,
    mode: Mode,
    write_back: WriteBack = WriteBack.NO,
    lock: Any = None,  # multiprocessing.Manager().Lock() is some crazy proxy
) -> bool:
    """Format file under `src` path. Return True if changed.

    If `write_back` is DIFF, write a diff to stdout. If it is YES, write reformatted
    code to the file.
    `mode` and `fast` options are passed to :func:`format_file_contents`.
    """
    if src.suffix == ".pyi":
        mode = replace(mode, is_pyi=True)

    then = datetime.utcfromtimestamp(src.stat().st_mtime)
    with open(src, "rb") as buf:
        src_contents, encoding, newline = decode_bytes(buf.read())
    try:
        dst_contents = format_file_contents(src_contents, fast=fast, mode=mode)
    except NothingChanged:
        return False

    if write_back == WriteBack.YES:
        with open(src, "w", encoding=encoding, newline=newline) as f:
            f.write(dst_contents)
    elif write_back in (WriteBack.DIFF, WriteBack.COLOR_DIFF):
        now = datetime.utcnow()
        src_name = f"{src}\t{then} +0000"
        dst_name = f"{src}\t{now} +0000"
        diff_contents = diff(src_contents, dst_contents, src_name, dst_name)

        if write_back == WriteBack.COLOR_DIFF:
            diff_contents = color_diff(diff_contents)

        with lock or nullcontext():
            f = io.TextIOWrapper(
                sys.stdout.buffer,
                encoding=encoding,
                newline=newline,
                write_through=True,
            )
            f = wrap_stream_for_windows(f)
            f.write(diff_contents)
            f.detach()

    return True


def format_stdin_to_stdout(
    fast: bool,
    *,
    content: Optional[str] = None,
    write_back: WriteBack = WriteBack.NO,
    mode: Mode,
) -> bool:
    """Format file on stdin. Return True if changed.

    If content is None, it's read from sys.stdin.

    If `write_back` is YES, write reformatted code back to stdout. If it is DIFF,
    write a diff to stdout. The `mode` argument is passed to
    :func:`format_file_contents`.
    """
    then = datetime.utcnow()

    if content is None:
        src, encoding, newline = decode_bytes(sys.stdin.buffer.read())
    else:
        src, encoding, newline = content, "utf-8", ""

    dst = src
    try:
        dst = format_file_contents(src, fast=fast, mode=mode)
        return True

    except NothingChanged:
        return False

    finally:
        f = io.TextIOWrapper(
            sys.stdout.buffer, encoding=encoding, newline=newline, write_through=True
        )
        if write_back == WriteBack.YES:
            # Make sure there's a newline after the content
            dst += "" if dst[-1] == "\n" else "\n"
            f.write(dst)
        elif write_back in (WriteBack.DIFF, WriteBack.COLOR_DIFF):
            now = datetime.utcnow()
            src_name = f"STDIN\t{then} +0000"
            dst_name = f"STDOUT\t{now} +0000"
            d = diff(src, dst, src_name, dst_name)
            if write_back == WriteBack.COLOR_DIFF:
                d = color_diff(d)
                f = wrap_stream_for_windows(f)
            f.write(d)
        f.detach()


def format_file_contents(src_contents: str, *, fast: bool, mode: Mode) -> FileContent:
    """Reformat contents of a file and return new contents.

    If `fast` is False, additionally confirm that the reformatted code is
    valid by calling :func:`assert_equivalent` and :func:`assert_stable` on it.
    `mode` is passed to :func:`format_str`.
    """
    if not src_contents.strip():
        raise NothingChanged

    dst_contents = format_str(src_contents, mode=mode)
    if src_contents == dst_contents:
        raise NothingChanged

    if not fast:
        assert_equivalent(src_contents, dst_contents)

        # Forced second pass to work around optional trailing commas (becoming
        # forced trailing commas on pass 2) interacting differently with optional
        # parentheses.  Admittedly ugly.
        dst_contents_pass2 = format_str(dst_contents, mode=mode)
        if dst_contents != dst_contents_pass2:
            dst_contents = dst_contents_pass2
            assert_equivalent(src_contents, dst_contents, pass_num=2)
            assert_stable(src_contents, dst_contents, mode=mode)
        # Note: no need to explicitly call `assert_stable` if `dst_contents` was
        # the same as `dst_contents_pass2`.
    return dst_contents


def format_str(src_contents: str, *, mode: Mode) -> FileContent:
    """Reformat a string and return new contents.

    `mode` determines formatting options, such as how many characters per line are
    allowed.  Example:

    >>> import black
    >>> print(black.format_str("def f(arg:str='')->None:...", mode=black.Mode()))
    def f(arg: str = "") -> None:
        ...

    A more complex example:

    >>> print(
    ...   black.format_str(
    ...     "def f(arg:str='')->None: hey",
    ...     mode=black.Mode(
    ...       target_versions={black.TargetVersion.PY36},
    ...       line_length=10,
    ...       string_normalization=False,
    ...       is_pyi=False,
    ...     ),
    ...   ),
    ... )
    def f(
        arg: str = '',
    ) -> None:
        hey

    """
    src_node = lib2to3_parse(src_contents.lstrip(), mode.target_versions)
    dst_contents = []
    future_imports = get_future_imports(src_node)
    if mode.target_versions:
        versions = mode.target_versions
    else:
        versions = detect_target_versions(src_node)
    normalize_fmt_off(src_node)
    lines = LineGenerator(
        mode=mode,
        remove_u_prefix="unicode_literals" in future_imports
        or supports_feature(versions, Feature.UNICODE_LITERALS),
    )
    elt = EmptyLineTracker(is_pyi=mode.is_pyi)
    empty_line = Line(mode=mode)
    after = 0
    split_line_features = {
        feature
        for feature in {Feature.TRAILING_COMMA_IN_CALL, Feature.TRAILING_COMMA_IN_DEF}
        if supports_feature(versions, feature)
    }
    for current_line in lines.visit(src_node):
        dst_contents.append(str(empty_line) * after)
        before, after = elt.maybe_empty_lines(current_line)
        dst_contents.append(str(empty_line) * before)
        for line in transform_line(
            current_line, mode=mode, features=split_line_features
        ):
            dst_contents.append(str(line))
    return "".join(dst_contents)


def decode_bytes(src: bytes) -> Tuple[FileContent, Encoding, NewLine]:
    """Return a tuple of (decoded_contents, encoding, newline).

    `newline` is either CRLF or LF but `decoded_contents` is decoded with
    universal newlines (i.e. only contains LF).
    """
    srcbuf = io.BytesIO(src)
    encoding, lines = tokenize.detect_encoding(srcbuf.readline)
    if not lines:
        return "", encoding, "\n"

    newline = "\r\n" if b"\r\n" == lines[0][-2:] else "\n"
    srcbuf.seek(0)
    with io.TextIOWrapper(srcbuf, encoding) as tiow:
        return tiow.read(), encoding, newline


def get_features_used(node: Node) -> Set[Feature]:
    """Return a set of (relatively) new Python features used in this file.

    Currently looking for:
    - f-strings;
    - underscores in numeric literals;
    - trailing commas after * or ** in function signatures and calls;
    - positional only arguments in function signatures and lambdas;
    - assignment expression;
    - relaxed decorator syntax;
    """
    features: Set[Feature] = set()
    for n in node.pre_order():
        if n.type == token.STRING:
            value_head = n.value[:2]  # type: ignore
            if value_head in {'f"', 'F"', "f'", "F'", "rf", "fr", "RF", "FR"}:
                features.add(Feature.F_STRINGS)

        elif n.type == token.NUMBER:
            if "_" in n.value:  # type: ignore
                features.add(Feature.NUMERIC_UNDERSCORES)

        elif n.type == token.SLASH:
            if n.parent and n.parent.type in {syms.typedargslist, syms.arglist}:
                features.add(Feature.POS_ONLY_ARGUMENTS)

        elif n.type == token.COLONEQUAL:
            features.add(Feature.ASSIGNMENT_EXPRESSIONS)

        elif n.type == syms.decorator:
            if len(n.children) > 1 and not is_simple_decorator_expression(
                n.children[1]
            ):
                features.add(Feature.RELAXED_DECORATORS)

        elif (
            n.type in {syms.typedargslist, syms.arglist}
            and n.children
            and n.children[-1].type == token.COMMA
        ):
            if n.type == syms.typedargslist:
                feature = Feature.TRAILING_COMMA_IN_DEF
            else:
                feature = Feature.TRAILING_COMMA_IN_CALL

            for ch in n.children:
                if ch.type in STARS:
                    features.add(feature)

                if ch.type == syms.argument:
                    for argch in ch.children:
                        if argch.type in STARS:
                            features.add(feature)

    return features


def detect_target_versions(node: Node) -> Set[TargetVersion]:
    """Detect the version to target based on the nodes used."""
    features = get_features_used(node)
    return {
        version for version in TargetVersion if features <= VERSION_TO_FEATURES[version]
    }


def get_future_imports(node: Node) -> Set[str]:
    """Return a set of __future__ imports in the file."""
    imports: Set[str] = set()

    def get_imports_from_children(children: List[LN]) -> Generator[str, None, None]:
        for child in children:
            if isinstance(child, Leaf):
                if child.type == token.NAME:
                    yield child.value

            elif child.type == syms.import_as_name:
                orig_name = child.children[0]
                assert isinstance(orig_name, Leaf), "Invalid syntax parsing imports"
                assert orig_name.type == token.NAME, "Invalid syntax parsing imports"
                yield orig_name.value

            elif child.type == syms.import_as_names:
                yield from get_imports_from_children(child.children)

            else:
                raise AssertionError("Invalid syntax parsing imports")

    for child in node.children:
        if child.type != syms.simple_stmt:
            break

        first_child = child.children[0]
        if isinstance(first_child, Leaf):
            # Continue looking if we see a docstring; otherwise stop.
            if (
                len(child.children) == 2
                and first_child.type == token.STRING
                and child.children[1].type == token.NEWLINE
            ):
                continue

            break

        elif first_child.type == syms.import_from:
            module_name = first_child.children[1]
            if not isinstance(module_name, Leaf) or module_name.value != "__future__":
                break

            imports |= set(get_imports_from_children(first_child.children[3:]))
        else:
            break

    return imports


def assert_equivalent(src: str, dst: str, *, pass_num: int = 1) -> None:
    """Raise AssertionError if `src` and `dst` aren't equivalent."""
    try:
        src_ast = parse_ast(src)
    except Exception as exc:
        raise AssertionError(
            "cannot use --safe with this file; failed to parse source file.  AST"
            f" error message: {exc}"
        )

    try:
        dst_ast = parse_ast(dst)
    except Exception as exc:
        log = dump_to_file("".join(traceback.format_tb(exc.__traceback__)), dst)
        raise AssertionError(
            f"INTERNAL ERROR: Black produced invalid code on pass {pass_num}: {exc}. "
            "Please report a bug on https://github.com/psf/black/issues.  "
            f"This invalid output might be helpful: {log}"
        ) from None

    src_ast_str = "\n".join(stringify_ast(src_ast))
    dst_ast_str = "\n".join(stringify_ast(dst_ast))
    if src_ast_str != dst_ast_str:
        log = dump_to_file(diff(src_ast_str, dst_ast_str, "src", "dst"))
        raise AssertionError(
            "INTERNAL ERROR: Black produced code that is not equivalent to the"
            f" source on pass {pass_num}.  Please report a bug on "
            f"https://github.com/psf/black/issues.  This diff might be helpful: {log}"
        ) from None


def assert_stable(src: str, dst: str, mode: Mode) -> None:
    """Raise AssertionError if `dst` reformats differently the second time."""
    newdst = format_str(dst, mode=mode)
    if dst != newdst:
        log = dump_to_file(
            str(mode),
            diff(src, dst, "source", "first pass"),
            diff(dst, newdst, "first pass", "second pass"),
        )
        raise AssertionError(
            "INTERNAL ERROR: Black produced different code on the second pass of the"
            " formatter.  Please report a bug on https://github.com/psf/black/issues."
            f"  This diff might be helpful: {log}"
        ) from None


@contextmanager
def nullcontext() -> Iterator[None]:
    """Return an empty context manager.

    To be used like `nullcontext` in Python 3.7.
    """
    yield


def patch_click() -> None:
    """Make Click not crash on Python 3.6 with LANG=C.

    On certain misconfigured environments, Python 3 selects the ASCII encoding as the
    default which restricts paths that it can access during the lifetime of the
    application.  Click refuses to work in this scenario by raising a RuntimeError.

    In case of Black the likelihood that non-ASCII characters are going to be used in
    file paths is minimal since it's Python source code.  Moreover, this crash was
    spurious on Python 3.7 thanks to PEP 538 and PEP 540.
    """
    try:
        from click import core
        from click import _unicodefun  # type: ignore
    except ModuleNotFoundError:
        return

    for module in (core, _unicodefun):
        if hasattr(module, "_verify_python3_env"):
            module._verify_python3_env = lambda: None  # type: ignore
        if hasattr(module, "_verify_python_env"):
            module._verify_python_env = lambda: None  # type: ignore


def patched_main() -> None:
    maybe_install_uvloop()
    freeze_support()
    patch_click()
    main()


if __name__ == "__main__":
    patched_main()
