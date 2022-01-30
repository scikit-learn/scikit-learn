"""Test cases for generating node-level dependencies (for fine-grained incremental checking)"""

import os
from collections import defaultdict

from typing import List, Tuple, Dict, Optional, Set
from typing_extensions import DefaultDict

from mypy import build, defaults
from mypy.modulefinder import BuildSource
from mypy.errors import CompileError
from mypy.nodes import MypyFile, Expression
from mypy.options import Options
from mypy.server.deps import get_dependencies
from mypy.test.config import test_temp_dir
from mypy.test.data import DataDrivenTestCase, DataSuite
from mypy.test.helpers import assert_string_arrays_equal, parse_options
from mypy.types import Type
from mypy.typestate import TypeState

# Only dependencies in these modules are dumped
dumped_modules = ['__main__', 'pkg', 'pkg.mod']


class GetDependenciesSuite(DataSuite):
    files = [
        'deps.test',
        'deps-types.test',
        'deps-generics.test',
        'deps-expressions.test',
        'deps-statements.test',
        'deps-classes.test',
    ]

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        src = '\n'.join(testcase.input)
        dump_all = '# __dump_all__' in src
        options = parse_options(src, testcase, incremental_step=1)
        if testcase.name.endswith('python2'):
            options.python_version = defaults.PYTHON2_VERSION
        options.use_builtins_fixtures = True
        options.show_traceback = True
        options.cache_dir = os.devnull
        options.export_types = True
        options.preserve_asts = True
        messages, files, type_map = self.build(src, options)
        a = messages
        if files is None or type_map is None:
            if not a:
                a = ['Unknown compile error (likely syntax error in test case or fixture)']
        else:
            deps: DefaultDict[str, Set[str]] = defaultdict(set)
            for module in files:
                if module in dumped_modules or dump_all and module not in ('abc',
                                                                           'typing',
                                                                           'mypy_extensions',
                                                                           'typing_extensions',
                                                                           'enum'):
                    new_deps = get_dependencies(files[module], type_map, options.python_version,
                                                options)
                    for source in new_deps:
                        deps[source].update(new_deps[source])

            TypeState.add_all_protocol_deps(deps)

            for source, targets in sorted(deps.items()):
                if source.startswith(('<enum', '<typing', '<mypy')):
                    # Remove noise.
                    continue
                line = '%s -> %s' % (source, ', '.join(sorted(targets)))
                # Clean up output a bit
                line = line.replace('__main__', 'm')
                a.append(line)

        assert_string_arrays_equal(
            testcase.output, a,
            'Invalid output ({}, line {})'.format(testcase.file,
                                                  testcase.line))

    def build(self,
              source: str,
              options: Options) -> Tuple[List[str],
                                         Optional[Dict[str, MypyFile]],
                                         Optional[Dict[Expression, Type]]]:
        try:
            result = build.build(sources=[BuildSource('main', None, source)],
                                 options=options,
                                 alt_lib_path=test_temp_dir)
        except CompileError as e:
            # TODO: Should perhaps not return None here.
            return e.messages, None, None
        return result.errors, result.files, result.types
