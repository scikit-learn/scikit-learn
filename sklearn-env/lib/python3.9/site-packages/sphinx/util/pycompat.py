"""
    sphinx.util.pycompat
    ~~~~~~~~~~~~~~~~~~~~

    Stuff for Python version compatibility.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import warnings
from typing import Any, Callable

from sphinx.deprecation import RemovedInSphinx60Warning

# ------------------------------------------------------------------------------
# Python 2/3 compatibility


# convert_with_2to3():
# support for running 2to3 over config files
def convert_with_2to3(filepath: str) -> str:
    warnings.warn('convert_with_2to3() is deprecated',
                  RemovedInSphinx60Warning, stacklevel=2)

    try:
        from lib2to3.pgen2.parse import ParseError
        from lib2to3.refactor import RefactoringTool, get_fixers_from_package
    except ImportError as exc:
        # python 3.9.0a6+ emits PendingDeprecationWarning for lib2to3.
        # Additionally, removal of the module is still discussed at PEP-594.
        # To support future python, this catches ImportError for lib2to3.
        raise SyntaxError from exc

    fixers = get_fixers_from_package('lib2to3.fixes')
    refactoring_tool = RefactoringTool(fixers)
    source = refactoring_tool._read_python_source(filepath)[0]
    try:
        tree = refactoring_tool.refactor_string(source, 'conf.py')
    except ParseError as err:
        # do not propagate lib2to3 exceptions
        lineno, offset = err.context[1]
        # try to match ParseError details with SyntaxError details

        raise SyntaxError(err.msg, (filepath, lineno, offset, err.value)) from err
    return str(tree)


def execfile_(filepath: str, _globals: Any, open: Callable = open) -> None:
    warnings.warn('execfile_() is deprecated',
                  RemovedInSphinx60Warning, stacklevel=2)
    from sphinx.util.osutil import fs_encoding
    with open(filepath, 'rb') as f:
        source = f.read()

    # compile to a code object, handle syntax errors
    filepath_enc = filepath.encode(fs_encoding)
    code = compile(source, filepath_enc, 'exec')
    exec(code, _globals)
