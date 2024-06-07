"""Build script for mypyc C runtime library unit tests.

The tests are written in C++ and use the Google Test framework.
"""

from __future__ import annotations

import os
import subprocess
import sys
from distutils.command.build_ext import build_ext
from distutils.core import Extension, setup
from typing import Any

kwargs: dict[str, Any]
if sys.platform == "darwin":
    kwargs = {"language": "c++"}
    compile_args = []
else:
    kwargs = {}
    compile_args = ["--std=c++11"]


class build_ext_custom(build_ext):
    def get_library_names(self):
        return ["gtest"]

    def run(self):
        gtest_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "external", "googletest")
        )

        os.makedirs(self.build_temp, exist_ok=True)

        # Build Google Test, the C++ framework we use for testing C code.
        # The source code for Google Test is copied to this repository.
        subprocess.check_call(
            [
                "make",
                "-f",
                os.path.join(gtest_dir, "make", "Makefile"),
                f"GTEST_DIR={gtest_dir}",
            ],
            cwd=self.build_temp,
        )

        self.library_dirs = [self.build_temp]

        return build_ext.run(self)


setup(
    name="test_capi",
    version="0.1",
    ext_modules=[
        Extension(
            "test_capi",
            [
                "test_capi.cc",
                "init.c",
                "int_ops.c",
                "float_ops.c",
                "list_ops.c",
                "exc_ops.c",
                "generic_ops.c",
            ],
            depends=["CPy.h", "mypyc_util.h", "pythonsupport.h"],
            extra_compile_args=["-Wno-unused-function", "-Wno-sign-compare"]
            + compile_args,
            libraries=["gtest"],
            include_dirs=["../external/googletest", "../external/googletest/include"],
            **kwargs,
        )
    ],
    cmdclass={"build_ext": build_ext_custom},
)
