# pylint: disable=missing-function-docstring
from __future__ import annotations

import argparse
import functools
import re
import sys
import sysconfig

from ._version import __version__
from .commands import get_cmake_dir, get_include, get_pkgconfig_dir

# This is the conditional used for os.path being posixpath
if "posix" in sys.builtin_module_names:
    from shlex import quote
elif "nt" in sys.builtin_module_names:
    # See https://github.com/mesonbuild/meson/blob/db22551ed9d2dd7889abea01cc1c7bba02bf1c75/mesonbuild/utils/universal.py#L1092-L1121
    # and the original documents:
    # https://docs.microsoft.com/en-us/cpp/c-language/parsing-c-command-line-arguments and
    # https://blogs.msdn.microsoft.com/twistylittlepassagesallalike/2011/04/23/everyone-quotes-command-line-arguments-the-wrong-way/
    UNSAFE = re.compile("[ \t\n\r]")

    def quote(s: str) -> str:
        if s and not UNSAFE.search(s):
            return s

        # Paths cannot contain a '"' on Windows, so we don't need to worry
        # about nuanced counting here.
        return f'"{s}\\"' if s.endswith("\\") else f'"{s}"'
else:

    def quote(s: str) -> str:
        return s


def print_includes() -> None:
    dirs = [
        sysconfig.get_path("include"),
        sysconfig.get_path("platinclude"),
        get_include(),
    ]

    # Make unique but preserve order
    unique_dirs = []
    for d in dirs:
        if d and d not in unique_dirs:
            unique_dirs.append(d)

    print(" ".join(quote(f"-I{d}") for d in unique_dirs))


def main() -> None:
    make_parser = functools.partial(argparse.ArgumentParser, allow_abbrev=False)
    if sys.version_info >= (3, 14):
        make_parser = functools.partial(make_parser, color=True, suggest_on_error=True)
    parser = make_parser()
    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
        help="Print the version and exit.",
    )
    parser.add_argument(
        "--includes",
        action="store_true",
        help="Include flags for both pybind11 and Python headers.",
    )
    parser.add_argument(
        "--cmakedir",
        action="store_true",
        help="Print the CMake module directory, ideal for setting -Dpybind11_ROOT in CMake.",
    )
    parser.add_argument(
        "--pkgconfigdir",
        action="store_true",
        help="Print the pkgconfig directory, ideal for setting $PKG_CONFIG_PATH.",
    )
    parser.add_argument(
        "--extension-suffix",
        action="store_true",
        help="Print the extension for a Python module",
    )
    args = parser.parse_args()
    if not sys.argv[1:]:
        parser.print_help()
    if args.includes:
        print_includes()
    if args.cmakedir:
        print(quote(get_cmake_dir()))
    if args.pkgconfigdir:
        print(quote(get_pkgconfig_dir()))
    if args.extension_suffix:
        print(sysconfig.get_config_var("EXT_SUFFIX"))


if __name__ == "__main__":
    main()
