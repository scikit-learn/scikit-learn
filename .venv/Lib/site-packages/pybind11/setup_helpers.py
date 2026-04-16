"""
This module provides helpers for C++11+ projects using pybind11.

LICENSE:

Copyright (c) 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>, All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# IMPORTANT: If you change this file in the pybind11 repo, also review
# setup_helpers.pyi for matching changes.
#
# If you copy this file in, you don't
# need the .pyi file; it's just an interface file for static type checkers.
from __future__ import annotations

import contextlib
import os
import platform
import shlex
import shutil
import sys
import sysconfig
import tempfile
import threading
import warnings
from functools import lru_cache
from pathlib import Path
from typing import (
    Any,
    Callable,
    Iterable,
    Iterator,
    List,
    Optional,
    Tuple,
    TypeVar,
    Union,
)

try:
    from setuptools import Extension as _Extension
    from setuptools.command.build_ext import build_ext as _build_ext
except ImportError:
    from distutils.command.build_ext import (  # type: ignore[assignment]
        build_ext as _build_ext,
    )
    from distutils.extension import Extension as _Extension  # type: ignore[assignment]

import distutils.ccompiler
import distutils.errors

WIN = sys.platform.startswith("win32") and "mingw" not in sysconfig.get_platform()
MACOS = sys.platform.startswith("darwin")
STD_TMPL = "/std:c++{}" if WIN else "-std=c++{}"


# It is recommended to use PEP 518 builds if using this module. However, this
# file explicitly supports being copied into a user's project directory
# standalone, and pulling pybind11 with the deprecated setup_requires feature.
# If you copy the file, remember to add it to your MANIFEST.in, and add the current
# directory into your path if it sits beside your setup.py.


class Pybind11Extension(_Extension):
    """
    Build a C++11+ Extension module with pybind11. This automatically adds the
    recommended flags when you init the extension and assumes C++ sources - you
    can further modify the options yourself.

    The customizations are:

    * ``/EHsc`` and ``/bigobj`` on Windows
    * ``stdlib=libc++`` on macOS
    * ``visibility=hidden`` and ``-g0`` on Unix

    Finally, you can set ``cxx_std`` via constructor or afterwards to enable
    flags for C++ std, and a few extra helper flags related to the C++ standard
    level. It is _highly_ recommended you either set this, or use the provided
    ``build_ext``, which will search for the highest supported extension for
    you if the ``cxx_std`` property is not set. Do not set the ``cxx_std``
    property more than once, as flags are added when you set it. Set the
    property to None to disable the addition of C++ standard flags.

    If you want to add pybind11 headers manually, for example for an exact
    git checkout, then set ``include_pybind11=False``.
    """

    # flags are prepended, so that they can be further overridden, e.g. by
    # ``extra_compile_args=["-g"]``.

    def _add_cflags(self, flags: list[str]) -> None:
        self.extra_compile_args[:0] = flags

    def _add_ldflags(self, flags: list[str]) -> None:
        self.extra_link_args[:0] = flags

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self._cxx_level = 0
        cxx_std = kwargs.pop("cxx_std", 0)

        if "language" not in kwargs:
            kwargs["language"] = "c++"

        include_pybind11 = kwargs.pop("include_pybind11", True)

        super().__init__(*args, **kwargs)

        # Include the installed package pybind11 headers
        if include_pybind11:
            # If using setup_requires, this fails the first time - that's okay
            try:
                import pybind11

                pyinc = pybind11.get_include()

                if pyinc not in self.include_dirs:
                    self.include_dirs.append(pyinc)
            except ModuleNotFoundError:
                pass

        self.cxx_std = cxx_std

        cflags = []
        if WIN:
            cflags += ["/EHsc", "/bigobj"]
        else:
            cflags += ["-fvisibility=hidden"]
            env_cflags = os.environ.get("CFLAGS", "")
            env_cppflags = os.environ.get("CPPFLAGS", "")
            c_cpp_flags = shlex.split(env_cflags) + shlex.split(env_cppflags)
            if not any(opt.startswith("-g") for opt in c_cpp_flags):
                cflags += ["-g0"]
        self._add_cflags(cflags)

    @property
    def cxx_std(self) -> int:
        """
        The CXX standard level. If set, will add the required flags. If left at
        0, it will trigger an automatic search when pybind11's build_ext is
        used. If None, will have no effect.  Besides just the flags, this may
        add a macos-min 10.9 or 10.14 flag if MACOSX_DEPLOYMENT_TARGET is
        unset.
        """
        return self._cxx_level

    @cxx_std.setter
    def cxx_std(self, level: int) -> None:
        if self._cxx_level:
            warnings.warn(
                "You cannot safely change the cxx_level after setting it!", stacklevel=2
            )

        # MSVC 2015 Update 3 and later only have 14 (and later 17) modes, so
        # force a valid flag here.
        if WIN and level == 11:
            level = 14

        self._cxx_level = level

        if not level:
            return

        cflags = [STD_TMPL.format(level)]
        ldflags = []

        if MACOS and "MACOSX_DEPLOYMENT_TARGET" not in os.environ:
            # C++17 requires a higher min version of macOS. An earlier version
            # (10.12 or 10.13) can be set manually via environment variable if
            # you are careful in your feature usage, but 10.14 is the safest
            # setting for general use. However, never set higher than the
            # current macOS version!
            current_macos = tuple(int(x) for x in platform.mac_ver()[0].split(".")[:2])
            desired_macos = (10, 9) if level < 17 else (10, 14)
            macos_string = ".".join(str(x) for x in min(current_macos, desired_macos))
            macosx_min = f"-mmacosx-version-min={macos_string}"
            cflags += [macosx_min]
            ldflags += [macosx_min]

        self._add_cflags(cflags)
        self._add_ldflags(ldflags)


# Just in case someone clever tries to multithread
tmp_chdir_lock = threading.Lock()


@contextlib.contextmanager
def tmp_chdir() -> Iterator[str]:
    "Prepare and enter a temporary directory, cleanup when done"

    # Threadsafe
    with tmp_chdir_lock:
        olddir = os.getcwd()
        try:
            tmpdir = tempfile.mkdtemp()
            os.chdir(tmpdir)
            yield tmpdir
        finally:
            os.chdir(olddir)
            shutil.rmtree(tmpdir)


# cf http://bugs.python.org/issue26689
def has_flag(compiler: Any, flag: str) -> bool:
    """
    Return the flag if a flag name is supported on the
    specified compiler, otherwise None (can be used as a boolean).
    If multiple flags are passed, return the first that matches.
    """

    with tmp_chdir():
        fname = Path("flagcheck.cpp")
        # Don't trigger -Wunused-parameter.
        fname.write_text("int main (int, char **) { return 0; }", encoding="utf-8")

        try:
            compiler.compile([str(fname)], extra_postargs=[flag])
        except distutils.errors.CompileError:
            return False
        return True


# Every call will cache the result
cpp_flag_cache = None


@lru_cache
def auto_cpp_level(compiler: Any) -> str | int:
    """
    Return the max supported C++ std level (17, 14, or 11). Returns latest on Windows.
    """

    if WIN:
        return "latest"

    levels = [17, 14, 11]

    for level in levels:
        if has_flag(compiler, STD_TMPL.format(level)):
            return level

    msg = "Unsupported compiler -- at least C++11 support is needed!"
    raise RuntimeError(msg)


class build_ext(_build_ext):  # noqa: N801
    """
    Customized build_ext that allows an auto-search for the highest supported
    C++ level for Pybind11Extension. This is only needed for the auto-search
    for now, and is completely optional otherwise.
    """

    def build_extensions(self) -> None:
        """
        Build extensions, injecting C++ std for Pybind11Extension if needed.
        """

        for ext in self.extensions:
            if hasattr(ext, "_cxx_level") and ext._cxx_level == 0:
                ext.cxx_std = auto_cpp_level(self.compiler)

        super().build_extensions()


def intree_extensions(
    paths: Iterable[str], package_dir: dict[str, str] | None = None
) -> list[Pybind11Extension]:
    """
    Generate Pybind11Extensions from source files directly located in a Python
    source tree.

    ``package_dir`` behaves as in ``setuptools.setup``.  If unset, the Python
    package root parent is determined as the first parent directory that does
    not contain an ``__init__.py`` file.
    """
    exts = []

    if package_dir is None:
        for path in paths:
            parent, _ = os.path.split(path)
            while os.path.exists(os.path.join(parent, "__init__.py")):
                parent, _ = os.path.split(parent)
            relname, _ = os.path.splitext(os.path.relpath(path, parent))
            qualified_name = relname.replace(os.path.sep, ".")
            exts.append(Pybind11Extension(qualified_name, [path]))
        return exts

    for path in paths:
        for prefix, parent in package_dir.items():
            if path.startswith(parent):
                relname, _ = os.path.splitext(os.path.relpath(path, parent))
                qualified_name = relname.replace(os.path.sep, ".")
                if prefix:
                    qualified_name = prefix + "." + qualified_name
                exts.append(Pybind11Extension(qualified_name, [path]))
                break
        else:
            msg = (
                f"path {path} is not a child of any of the directories listed "
                f"in 'package_dir' ({package_dir})"
            )
            raise ValueError(msg)

    return exts


def naive_recompile(obj: str, src: str) -> bool:
    """
    This will recompile only if the source file changes. It does not check
    header files, so a more advanced function or Ccache is better if you have
    editable header files in your package.
    """
    return os.stat(obj).st_mtime < os.stat(src).st_mtime


def no_recompile(obg: str, src: str) -> bool:  # noqa: ARG001
    """
    This is the safest but slowest choice (and is the default) - will always
    recompile sources.
    """
    return True


S = TypeVar("S", bound="ParallelCompile")

CCompilerMethod = Callable[
    [
        distutils.ccompiler.CCompiler,
        List[str],
        Optional[str],
        Optional[List[Union[Tuple[str], Tuple[str, Optional[str]]]]],
        Optional[List[str]],
        bool,
        Optional[List[str]],
        Optional[List[str]],
        Optional[List[str]],
    ],
    List[str],
]


# Optional parallel compile utility
# inspired by: http://stackoverflow.com/questions/11013851/speeding-up-build-process-with-distutils
# and: https://github.com/tbenthompson/cppimport/blob/stable/cppimport/build_module.py
# and NumPy's parallel distutils module:
#              https://github.com/numpy/numpy/blob/master/numpy/distutils/ccompiler.py
class ParallelCompile:
    """
    Make a parallel compile function. Inspired by
    numpy.distutils.ccompiler.CCompiler.compile and cppimport.

    This takes several arguments that allow you to customize the compile
    function created:

    envvar:
        Set an environment variable to control the compilation threads, like
        NPY_NUM_BUILD_JOBS
    default:
        0 will automatically multithread, or 1 will only multithread if the
        envvar is set.
    max:
        The limit for automatic multithreading if non-zero
    needs_recompile:
        A function of (obj, src) that returns True when recompile is needed.  No
        effect in isolated mode; use ccache instead, see
        https://github.com/matplotlib/matplotlib/issues/1507/

    To use::

        ParallelCompile("NPY_NUM_BUILD_JOBS").install()

    or::

        with ParallelCompile("NPY_NUM_BUILD_JOBS"):
            setup(...)

    By default, this assumes all files need to be recompiled. A smarter
    function can be provided via needs_recompile.  If the output has not yet
    been generated, the compile will always run, and this function is not
    called.
    """

    __slots__ = ("envvar", "default", "max", "_old", "needs_recompile")

    def __init__(
        self,
        envvar: str | None = None,
        default: int = 0,
        max: int = 0,  # pylint: disable=redefined-builtin
        needs_recompile: Callable[[str, str], bool] = no_recompile,
    ) -> None:
        self.envvar = envvar
        self.default = default
        self.max = max
        self.needs_recompile = needs_recompile
        self._old: list[CCompilerMethod] = []

    def function(self) -> CCompilerMethod:
        """
        Builds a function object usable as distutils.ccompiler.CCompiler.compile.
        """

        def compile_function(
            compiler: distutils.ccompiler.CCompiler,
            sources: list[str],
            output_dir: str | None = None,
            macros: list[tuple[str] | tuple[str, str | None]] | None = None,
            include_dirs: list[str] | None = None,
            debug: bool = False,
            extra_preargs: list[str] | None = None,
            extra_postargs: list[str] | None = None,
            depends: list[str] | None = None,
        ) -> Any:
            # These lines are directly from distutils.ccompiler.CCompiler
            macros, objects, extra_postargs, pp_opts, build = compiler._setup_compile(  # type: ignore[attr-defined]
                output_dir, macros, include_dirs, sources, depends, extra_postargs
            )
            cc_args = compiler._get_cc_args(pp_opts, debug, extra_preargs)  # type: ignore[attr-defined]

            # The number of threads; start with default.
            threads = self.default

            # Determine the number of compilation threads, unless set by an environment variable.
            if self.envvar is not None:
                threads = int(os.environ.get(self.envvar, self.default))

            def _single_compile(obj: Any) -> None:
                try:
                    src, ext = build[obj]
                except KeyError:
                    return

                if not os.path.exists(obj) or self.needs_recompile(obj, src):
                    compiler._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)  # type: ignore[attr-defined]

            try:
                # Importing .synchronize checks for platforms that have some multiprocessing
                # capabilities but lack semaphores, such as AWS Lambda and Android Termux.
                import multiprocessing.synchronize
                from multiprocessing.pool import ThreadPool
            except ImportError:
                threads = 1

            if threads == 0:
                try:
                    threads = multiprocessing.cpu_count()
                    threads = self.max if self.max and self.max < threads else threads
                except NotImplementedError:
                    threads = 1

            if threads > 1:
                with ThreadPool(threads) as pool:
                    for _ in pool.imap_unordered(_single_compile, objects):
                        pass
            else:
                for ob in objects:
                    _single_compile(ob)

            return objects

        return compile_function

    def install(self: S) -> S:
        """
        Installs the compile function into distutils.ccompiler.CCompiler.compile.
        """
        distutils.ccompiler.CCompiler.compile = self.function()  # type: ignore[assignment]
        return self

    def __enter__(self: S) -> S:
        self._old.append(distutils.ccompiler.CCompiler.compile)
        return self.install()

    def __exit__(self, *args: Any) -> None:
        distutils.ccompiler.CCompiler.compile = self._old.pop()  # type: ignore[assignment]
