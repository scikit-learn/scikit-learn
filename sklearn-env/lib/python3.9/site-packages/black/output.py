"""Nice output for Black.

The double calls are for patching purposes in tests.
"""

from typing import Any, Optional
from mypy_extensions import mypyc_attr
import tempfile

from click import echo, style


def _out(message: Optional[str] = None, nl: bool = True, **styles: Any) -> None:
    if message is not None:
        if "bold" not in styles:
            styles["bold"] = True
        message = style(message, **styles)
    echo(message, nl=nl, err=True)


def _err(message: Optional[str] = None, nl: bool = True, **styles: Any) -> None:
    if message is not None:
        if "fg" not in styles:
            styles["fg"] = "red"
        message = style(message, **styles)
    echo(message, nl=nl, err=True)


def out(message: Optional[str] = None, nl: bool = True, **styles: Any) -> None:
    _out(message, nl=nl, **styles)


def err(message: Optional[str] = None, nl: bool = True, **styles: Any) -> None:
    _err(message, nl=nl, **styles)


def diff(a: str, b: str, a_name: str, b_name: str) -> str:
    """Return a unified diff string between strings `a` and `b`."""
    import difflib

    a_lines = [line for line in a.splitlines(keepends=True)]
    b_lines = [line for line in b.splitlines(keepends=True)]
    diff_lines = []
    for line in difflib.unified_diff(
        a_lines, b_lines, fromfile=a_name, tofile=b_name, n=5
    ):
        # Work around https://bugs.python.org/issue2142
        # See:
        # https://www.gnu.org/software/diffutils/manual/html_node/Incomplete-Lines.html
        if line[-1] == "\n":
            diff_lines.append(line)
        else:
            diff_lines.append(line + "\n")
            diff_lines.append("\\ No newline at end of file\n")
    return "".join(diff_lines)


def color_diff(contents: str) -> str:
    """Inject the ANSI color codes to the diff."""
    lines = contents.split("\n")
    for i, line in enumerate(lines):
        if line.startswith("+++") or line.startswith("---"):
            line = "\033[1;37m" + line + "\033[0m"  # bold white, reset
        elif line.startswith("@@"):
            line = "\033[36m" + line + "\033[0m"  # cyan, reset
        elif line.startswith("+"):
            line = "\033[32m" + line + "\033[0m"  # green, reset
        elif line.startswith("-"):
            line = "\033[31m" + line + "\033[0m"  # red, reset
        lines[i] = line
    return "\n".join(lines)


@mypyc_attr(patchable=True)
def dump_to_file(*output: str, ensure_final_newline: bool = True) -> str:
    """Dump `output` to a temporary file. Return path to the file."""
    with tempfile.NamedTemporaryFile(
        mode="w", prefix="blk_", suffix=".log", delete=False, encoding="utf8"
    ) as f:
        for lines in output:
            f.write(lines)
            if ensure_final_newline and lines and lines[-1] != "\n":
                f.write("\n")
    return f.name
