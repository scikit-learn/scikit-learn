"""
pygments.lexers.agile
~~~~~~~~~~~~~~~~~~~~~

Just export lexer classes previously contained in this module.

:copyright: Copyright 2006-2025 by the Pygments team, see AUTHORS.
:license: BSD, see LICENSE for details.
"""

# ruff: noqa: F401

from pygments.lexers.d import CrocLexer, MiniDLexer
from pygments.lexers.factor import FactorLexer
from pygments.lexers.iolang import IoLexer
from pygments.lexers.jvm import ClojureLexer, IokeLexer
from pygments.lexers.lisp import SchemeLexer
from pygments.lexers.perl import Perl6Lexer, PerlLexer
from pygments.lexers.python import (
    DgLexer,
    Python3Lexer,
    Python3TracebackLexer,
    PythonConsoleLexer,
    PythonLexer,
    PythonTracebackLexer,
)
from pygments.lexers.ruby import FancyLexer, RubyConsoleLexer, RubyLexer
from pygments.lexers.scripting import LuaLexer, MoonScriptLexer
from pygments.lexers.tcl import TclLexer

__all__ = []
