"""Test cases for the type checker: exporting inferred types"""

import re

from mypy import build
from mypy.modulefinder import BuildSource
from mypy.test.config import test_temp_dir
from mypy.test.data import DataDrivenTestCase, DataSuite
from mypy.test.helpers import assert_string_arrays_equal
from mypy.test.visitors import SkippedNodeSearcher, ignore_node
from mypy.util import short_type
from mypy.nodes import NameExpr
from mypy.errors import CompileError
from mypy.options import Options


class TypeExportSuite(DataSuite):
    required_out_section = True
    files = ['typexport-basic.test']

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        try:
            line = testcase.input[0]
            mask = ''
            if line.startswith('##'):
                mask = '(' + line[2:].strip() + ')$'

            src = '\n'.join(testcase.input)
            options = Options()
            options.strict_optional = False  # TODO: Enable strict optional checking
            options.use_builtins_fixtures = True
            options.show_traceback = True
            options.export_types = True
            options.preserve_asts = True
            result = build.build(sources=[BuildSource('main', None, src)],
                                 options=options,
                                 alt_lib_path=test_temp_dir)
            a = result.errors
            map = result.types
            nodes = map.keys()

            # Ignore NameExpr nodes of variables with explicit (trivial) types
            # to simplify output.
            searcher = SkippedNodeSearcher()
            for file in result.files.values():
                file.accept(searcher)
            ignored = searcher.nodes

            # Filter nodes that should be included in the output.
            keys = []
            for node in nodes:
                if node.line is not None and node.line != -1 and map[node]:
                    if ignore_node(node) or node in ignored:
                        continue
                    if (re.match(mask, short_type(node))
                            or (isinstance(node, NameExpr)
                                and re.match(mask, node.name))):
                        # Include node in output.
                        keys.append(node)

            for key in sorted(keys,
                              key=lambda n: (n.line, short_type(n),
                                             str(n) + str(map[n]))):
                ts = str(map[key]).replace('*', '')  # Remove erased tags
                ts = ts.replace('__main__.', '')
                a.append('{}({}) : {}'.format(short_type(key), key.line, ts))
        except CompileError as e:
            a = e.messages
        assert_string_arrays_equal(
            testcase.output, a,
            'Invalid type checker output ({}, line {})'.format(testcase.file,
                                                               testcase.line))
