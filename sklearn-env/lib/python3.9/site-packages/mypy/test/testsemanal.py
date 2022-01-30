"""Semantic analyzer test cases"""

import os.path

from typing import Dict, List

from mypy import build
from mypy.modulefinder import BuildSource
from mypy.defaults import PYTHON3_VERSION
from mypy.test.helpers import (
    assert_string_arrays_equal, normalize_error_messages, testfile_pyversion, parse_options
)
from mypy.test.data import DataDrivenTestCase, DataSuite
from mypy.test.config import test_temp_dir
from mypy.errors import CompileError
from mypy.nodes import TypeInfo
from mypy.options import Options


# Semantic analyzer test cases: dump parse tree

# Semantic analysis test case description files.
semanal_files = [
    'semanal-basic.test',
    'semanal-expressions.test',
    'semanal-classes.test',
    'semanal-types.test',
    'semanal-typealiases.test',
    'semanal-modules.test',
    'semanal-statements.test',
    'semanal-abstractclasses.test',
    'semanal-namedtuple.test',
    'semanal-typeddict.test',
    'semenal-literal.test',
    'semanal-classvar.test',
    'semanal-python2.test',
    'semanal-lambda.test',
]


def get_semanal_options(program_text: str, testcase: DataDrivenTestCase) -> Options:
    options = parse_options(program_text, testcase, 1)
    options.use_builtins_fixtures = True
    options.semantic_analysis_only = True
    options.show_traceback = True
    options.python_version = PYTHON3_VERSION
    return options


class SemAnalSuite(DataSuite):
    files = semanal_files
    native_sep = True

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        test_semanal(testcase)


def test_semanal(testcase: DataDrivenTestCase) -> None:
    """Perform a semantic analysis test case.

    The testcase argument contains a description of the test case
    (inputs and output).
    """

    try:
        src = '\n'.join(testcase.input)
        options = get_semanal_options(src, testcase)
        options.python_version = testfile_pyversion(testcase.file)
        result = build.build(sources=[BuildSource('main', None, src)],
                             options=options,
                             alt_lib_path=test_temp_dir)
        a = result.errors
        if a:
            raise CompileError(a)
        # Include string representations of the source files in the actual
        # output.
        for fnam in sorted(result.files.keys()):
            f = result.files[fnam]
            # Omit the builtins module and files with a special marker in the
            # path.
            # TODO the test is not reliable
            if (not f.path.endswith((os.sep + 'builtins.pyi',
                                     'typing.pyi',
                                     'mypy_extensions.pyi',
                                     'typing_extensions.pyi',
                                     'abc.pyi',
                                     'collections.pyi',
                                     'sys.pyi'))
                    and not os.path.basename(f.path).startswith('_')
                    and not os.path.splitext(
                        os.path.basename(f.path))[0].endswith('_')):
                a += str(f).split('\n')
    except CompileError as e:
        a = e.messages
    if testcase.normalize_output:
        a = normalize_error_messages(a)
    assert_string_arrays_equal(
        testcase.output, a,
        'Invalid semantic analyzer output ({}, line {})'.format(testcase.file,
                                                                testcase.line))


# Semantic analyzer error test cases

class SemAnalErrorSuite(DataSuite):
    files = ['semanal-errors.test']

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        test_semanal_error(testcase)


def test_semanal_error(testcase: DataDrivenTestCase) -> None:
    """Perform a test case."""

    try:
        src = '\n'.join(testcase.input)
        res = build.build(sources=[BuildSource('main', None, src)],
                          options=get_semanal_options(src, testcase),
                          alt_lib_path=test_temp_dir)
        a = res.errors
        assert a, 'No errors reported in {}, line {}'.format(testcase.file, testcase.line)
    except CompileError as e:
        # Verify that there was a compile error and that the error messages
        # are equivalent.
        a = e.messages
    if testcase.normalize_output:
        a = normalize_error_messages(a)
    assert_string_arrays_equal(
        testcase.output, a,
        'Invalid compiler output ({}, line {})'.format(testcase.file, testcase.line))


# SymbolNode table export test cases

class SemAnalSymtableSuite(DataSuite):
    required_out_section = True
    files = ['semanal-symtable.test']

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        """Perform a test case."""
        try:
            # Build test case input.
            src = '\n'.join(testcase.input)
            result = build.build(sources=[BuildSource('main', None, src)],
                                 options=get_semanal_options(src, testcase),
                                 alt_lib_path=test_temp_dir)
            # The output is the symbol table converted into a string.
            a = result.errors
            if a:
                raise CompileError(a)
            for f in sorted(result.files.keys()):
                if f not in ('builtins', 'typing', 'abc'):
                    a.append('{}:'.format(f))
                    for s in str(result.files[f].names).split('\n'):
                        a.append('  ' + s)
        except CompileError as e:
            a = e.messages
        assert_string_arrays_equal(
            testcase.output, a,
            'Invalid semantic analyzer output ({}, line {})'.format(
                testcase.file, testcase.line))


# Type info export test cases
class SemAnalTypeInfoSuite(DataSuite):
    required_out_section = True
    files = ['semanal-typeinfo.test']

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        """Perform a test case."""
        try:
            # Build test case input.
            src = '\n'.join(testcase.input)
            result = build.build(sources=[BuildSource('main', None, src)],
                                 options=get_semanal_options(src, testcase),
                                 alt_lib_path=test_temp_dir)
            a = result.errors
            if a:
                raise CompileError(a)

            # Collect all TypeInfos in top-level modules.
            typeinfos = TypeInfoMap()
            for f in result.files.values():
                for n in f.names.values():
                    if isinstance(n.node, TypeInfo):
                        assert n.fullname is not None
                        typeinfos[n.fullname] = n.node

            # The output is the symbol table converted into a string.
            a = str(typeinfos).split('\n')
        except CompileError as e:
            a = e.messages
        assert_string_arrays_equal(
            testcase.output, a,
            'Invalid semantic analyzer output ({}, line {})'.format(
                testcase.file, testcase.line))


class TypeInfoMap(Dict[str, TypeInfo]):
    def __str__(self) -> str:
        a: List[str] = ["TypeInfoMap("]
        for x, y in sorted(self.items()):
            if isinstance(x, str) and (not x.startswith('builtins.') and
                                       not x.startswith('typing.') and
                                       not x.startswith('abc.')):
                ti = ('\n' + '  ').join(str(y).split('\n'))
                a.append('  {} : {}'.format(x, ti))
        a[-1] += ')'
        return '\n'.join(a)
