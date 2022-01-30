"""Type checker test cases"""

import os
import re
import sys

from typing import Dict, List, Set, Tuple

from mypy import build
from mypy.build import Graph
from mypy.modulefinder import BuildSource, SearchPaths, FindModuleCache
from mypy.test.config import test_temp_dir, test_data_prefix
from mypy.test.data import (
    DataDrivenTestCase, DataSuite, FileOperation, module_from_path
)
from mypy.test.helpers import (
    assert_string_arrays_equal, normalize_error_messages, assert_module_equivalence,
    update_testcase_output, parse_options,
    assert_target_equivalence, check_test_output_files, perform_file_operations,
)
from mypy.errors import CompileError
from mypy.semanal_main import core_modules


# List of files that contain test case descriptions.
typecheck_files = [
    'check-basic.test',
    'check-union-or-syntax.test',
    'check-callable.test',
    'check-classes.test',
    'check-statements.test',
    'check-generics.test',
    'check-dynamic-typing.test',
    'check-inference.test',
    'check-inference-context.test',
    'check-kwargs.test',
    'check-overloading.test',
    'check-type-checks.test',
    'check-abstract.test',
    'check-multiple-inheritance.test',
    'check-super.test',
    'check-modules.test',
    'check-typevar-values.test',
    'check-unsupported.test',
    'check-unreachable-code.test',
    'check-unions.test',
    'check-isinstance.test',
    'check-lists.test',
    'check-namedtuple.test',
    'check-narrowing.test',
    'check-typeddict.test',
    'check-type-aliases.test',
    'check-ignore.test',
    'check-type-promotion.test',
    'check-semanal-error.test',
    'check-flags.test',
    'check-incremental.test',
    'check-serialize.test',
    'check-bound.test',
    'check-optional.test',
    'check-fastparse.test',
    'check-warnings.test',
    'check-async-await.test',
    'check-newtype.test',
    'check-class-namedtuple.test',
    'check-selftype.test',
    'check-python2.test',
    'check-columns.test',
    'check-functions.test',
    'check-tuples.test',
    'check-expressions.test',
    'check-generic-subtyping.test',
    'check-varargs.test',
    'check-newsyntax.test',
    'check-protocols.test',
    'check-underscores.test',
    'check-classvar.test',
    'check-enum.test',
    'check-incomplete-fixture.test',
    'check-custom-plugin.test',
    'check-default-plugin.test',
    'check-attr.test',
    'check-ctypes.test',
    'check-dataclasses.test',
    'check-final.test',
    'check-redefine.test',
    'check-literal.test',
    'check-newsemanal.test',
    'check-inline-config.test',
    'check-reports.test',
    'check-errorcodes.test',
    'check-annotated.test',
    'check-parameter-specification.test',
    'check-generic-alias.test',
    'check-typeguard.test',
    'check-functools.test',
    'check-singledispatch.test',
    'check-slots.test',
    'check-formatting.test',
]

# Tests that use Python 3.8-only AST features (like expression-scoped ignores):
if sys.version_info >= (3, 8):
    typecheck_files.append('check-python38.test')
if sys.version_info >= (3, 9):
    typecheck_files.append('check-python39.test')
if sys.version_info >= (3, 10):
    typecheck_files.append('check-python310.test')

# Special tests for platforms with case-insensitive filesystems.
if sys.platform in ('darwin', 'win32'):
    typecheck_files.extend(['check-modules-case.test'])


class TypeCheckSuite(DataSuite):
    files = typecheck_files

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        incremental = ('incremental' in testcase.name.lower()
                       or 'incremental' in testcase.file
                       or 'serialize' in testcase.file)
        if incremental:
            # Incremental tests are run once with a cold cache, once with a warm cache.
            # Expect success on first run, errors from testcase.output (if any) on second run.
            num_steps = max([2] + list(testcase.output2.keys()))
            # Check that there are no file changes beyond the last run (they would be ignored).
            for dn, dirs, files in os.walk(os.curdir):
                for file in files:
                    m = re.search(r'\.([2-9])$', file)
                    if m and int(m.group(1)) > num_steps:
                        raise ValueError(
                            'Output file {} exists though test case only has {} runs'.format(
                                file, num_steps))
            steps = testcase.find_steps()
            for step in range(1, num_steps + 1):
                idx = step - 2
                ops = steps[idx] if idx < len(steps) and idx >= 0 else []
                self.run_case_once(testcase, ops, step)
        else:
            self.run_case_once(testcase)

    def run_case_once(self, testcase: DataDrivenTestCase,
                      operations: List[FileOperation] = [],
                      incremental_step: int = 0) -> None:
        original_program_text = '\n'.join(testcase.input)
        module_data = self.parse_module(original_program_text, incremental_step)

        # Unload already loaded plugins, they may be updated.
        for file, _ in testcase.files:
            module = module_from_path(file)
            if module.endswith('_plugin') and module in sys.modules:
                del sys.modules[module]
        if incremental_step == 0 or incremental_step == 1:
            # In run 1, copy program text to program file.
            for module_name, program_path, program_text in module_data:
                if module_name == '__main__':
                    with open(program_path, 'w', encoding='utf8') as f:
                        f.write(program_text)
                    break
        elif incremental_step > 1:
            # In runs 2+, copy *.[num] files to * files.
            perform_file_operations(operations)

        # Parse options after moving files (in case mypy.ini is being moved).
        options = parse_options(original_program_text, testcase, incremental_step)
        options.use_builtins_fixtures = True
        options.show_traceback = True

        # Enable some options automatically based on test file name.
        if 'optional' in testcase.file:
            options.strict_optional = True
        if 'columns' in testcase.file:
            options.show_column_numbers = True
        if 'errorcodes' in testcase.file:
            options.show_error_codes = True

        if incremental_step and options.incremental:
            # Don't overwrite # flags: --no-incremental in incremental test cases
            options.incremental = True
        else:
            options.incremental = False
            # Don't waste time writing cache unless we are specifically looking for it
            if not testcase.writescache:
                options.cache_dir = os.devnull

        sources = []
        for module_name, program_path, program_text in module_data:
            # Always set to none so we're forced to reread the module in incremental mode
            sources.append(BuildSource(program_path, module_name,
                                       None if incremental_step else program_text))

        plugin_dir = os.path.join(test_data_prefix, 'plugins')
        sys.path.insert(0, plugin_dir)

        res = None
        try:
            res = build.build(sources=sources,
                              options=options,
                              alt_lib_path=test_temp_dir)
            a = res.errors
        except CompileError as e:
            a = e.messages
        finally:
            assert sys.path[0] == plugin_dir
            del sys.path[0]

        if testcase.normalize_output:
            a = normalize_error_messages(a)

        # Make sure error messages match
        if incremental_step == 0:
            # Not incremental
            msg = 'Unexpected type checker output ({}, line {})'
            output = testcase.output
        elif incremental_step == 1:
            msg = 'Unexpected type checker output in incremental, run 1 ({}, line {})'
            output = testcase.output
        elif incremental_step > 1:
            msg = ('Unexpected type checker output in incremental, run {}'.format(
                incremental_step) + ' ({}, line {})')
            output = testcase.output2.get(incremental_step, [])
        else:
            raise AssertionError()

        if output != a and testcase.config.getoption('--update-data', False):
            update_testcase_output(testcase, a)
        assert_string_arrays_equal(output, a, msg.format(testcase.file, testcase.line))

        if res:
            if options.cache_dir != os.devnull:
                self.verify_cache(module_data, res.errors, res.manager, res.graph)

            name = 'targets'
            if incremental_step:
                name += str(incremental_step + 1)
            expected = testcase.expected_fine_grained_targets.get(incremental_step + 1)
            actual = res.manager.processed_targets
            # Skip the initial builtin cycle.
            actual = [t for t in actual
                      if not any(t.startswith(mod)
                                 for mod in core_modules + ['mypy_extensions'])]
            if expected is not None:
                assert_target_equivalence(name, expected, actual)
            if incremental_step > 1:
                suffix = '' if incremental_step == 2 else str(incremental_step - 1)
                expected_rechecked = testcase.expected_rechecked_modules.get(incremental_step - 1)
                if expected_rechecked is not None:
                    assert_module_equivalence(
                        'rechecked' + suffix,
                        expected_rechecked, res.manager.rechecked_modules)
                expected_stale = testcase.expected_stale_modules.get(incremental_step - 1)
                if expected_stale is not None:
                    assert_module_equivalence(
                        'stale' + suffix,
                        expected_stale, res.manager.stale_modules)

        if testcase.output_files:
            check_test_output_files(testcase, incremental_step, strip_prefix='tmp/')

    def verify_cache(self, module_data: List[Tuple[str, str, str]], a: List[str],
                     manager: build.BuildManager, graph: Graph) -> None:
        # There should be valid cache metadata for each module except
        # for those that had an error in themselves or one of their
        # dependencies.
        error_paths = self.find_error_message_paths(a)
        busted_paths = {m.path for id, m in manager.modules.items()
                        if graph[id].transitive_error}
        modules = self.find_module_files(manager)
        modules.update({module_name: path for module_name, path, text in module_data})
        missing_paths = self.find_missing_cache_files(modules, manager)
        # We would like to assert error_paths.issubset(busted_paths)
        # but this runs into trouble because while some 'notes' are
        # really errors that cause an error to be marked, many are
        # just notes attached to other errors.
        assert error_paths or not busted_paths, "Some modules reported error despite no errors"
        if not missing_paths == busted_paths:
            raise AssertionError("cache data discrepancy %s != %s" %
                                 (missing_paths, busted_paths))
        assert os.path.isfile(os.path.join(manager.options.cache_dir, ".gitignore"))
        cachedir_tag = os.path.join(manager.options.cache_dir, "CACHEDIR.TAG")
        assert os.path.isfile(cachedir_tag)
        with open(cachedir_tag) as f:
            assert f.read().startswith("Signature: 8a477f597d28d172789f06886806bc55")

    def find_error_message_paths(self, a: List[str]) -> Set[str]:
        hits = set()
        for line in a:
            m = re.match(r'([^\s:]+):(\d+:)?(\d+:)? (error|warning|note):', line)
            if m:
                p = m.group(1)
                hits.add(p)
        return hits

    def find_module_files(self, manager: build.BuildManager) -> Dict[str, str]:
        modules = {}
        for id, module in manager.modules.items():
            modules[id] = module.path
        return modules

    def find_missing_cache_files(self, modules: Dict[str, str],
                                 manager: build.BuildManager) -> Set[str]:
        ignore_errors = True
        missing = {}
        for id, path in modules.items():
            meta = build.find_cache_meta(id, path, manager)
            if not build.validate_meta(meta, id, path, ignore_errors, manager):
                missing[id] = path
        return set(missing.values())

    def parse_module(self,
                     program_text: str,
                     incremental_step: int = 0) -> List[Tuple[str, str, str]]:
        """Return the module and program names for a test case.

        Normally, the unit tests will parse the default ('__main__')
        module and follow all the imports listed there. You can override
        this behavior and instruct the tests to check multiple modules
        by using a comment like this in the test case input:

          # cmd: mypy -m foo.bar foo.baz

        You can also use `# cmdN:` to have a different cmd for incremental
        step N (2, 3, ...).

        Return a list of tuples (module name, file name, program text).
        """
        m = re.search('# cmd: mypy -m ([a-zA-Z0-9_. ]+)$', program_text, flags=re.MULTILINE)
        if incremental_step > 1:
            alt_regex = '# cmd{}: mypy -m ([a-zA-Z0-9_. ]+)$'.format(incremental_step)
            alt_m = re.search(alt_regex, program_text, flags=re.MULTILINE)
            if alt_m is not None:
                # Optionally return a different command if in a later step
                # of incremental mode, otherwise default to reusing the
                # original cmd.
                m = alt_m

        if m:
            # The test case wants to use a non-default main
            # module. Look up the module and give it as the thing to
            # analyze.
            module_names = m.group(1)
            out = []
            search_paths = SearchPaths((test_temp_dir,), (), (), ())
            cache = FindModuleCache(search_paths, fscache=None, options=None)
            for module_name in module_names.split(' '):
                path = cache.find_module(module_name)
                assert isinstance(path, str), "Can't find ad hoc case file: %s" % module_name
                with open(path, encoding='utf8') as f:
                    program_text = f.read()
                out.append((module_name, path, program_text))
            return out
        else:
            return [('__main__', 'main', program_text)]
