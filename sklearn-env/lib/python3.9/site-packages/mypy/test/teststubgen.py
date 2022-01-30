import io
import os.path
import shutil
import sys
import tempfile
import re
import unittest
from types import ModuleType

from typing import Any, List, Tuple, Optional

from mypy.test.helpers import (
    assert_equal, assert_string_arrays_equal, local_sys_path_set
)
from mypy.test.data import DataSuite, DataDrivenTestCase
from mypy.errors import CompileError
from mypy.stubgen import (
    generate_stubs, parse_options, Options, collect_build_targets,
    mypy_options, is_blacklisted_path, is_non_library_module
)
from mypy.stubutil import walk_packages, remove_misplaced_type_comments, common_dir_prefix
from mypy.stubgenc import (
    generate_c_type_stub, infer_method_sig, generate_c_function_stub, generate_c_property_stub
)
from mypy.stubdoc import (
    parse_signature, parse_all_signatures, build_signature, find_unique_signatures,
    infer_sig_from_docstring, infer_prop_type_from_docstring, FunctionSig, ArgSig,
    infer_arg_sig_from_anon_docstring, is_valid_type
)
from mypy.moduleinspect import ModuleInspect, InspectError


class StubgenCmdLineSuite(unittest.TestCase):
    """Test cases for processing command-line options and finding files."""

    @unittest.skipIf(sys.platform == 'win32', "clean up fails on Windows")
    def test_files_found(self) -> None:
        current = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            try:
                os.chdir(tmp)
                os.mkdir('subdir')
                self.make_file('subdir', 'a.py')
                self.make_file('subdir', 'b.py')
                os.mkdir(os.path.join('subdir', 'pack'))
                self.make_file('subdir', 'pack', '__init__.py')
                opts = parse_options(['subdir'])
                py_mods, c_mods = collect_build_targets(opts, mypy_options(opts))
                assert_equal(c_mods, [])
                files = {mod.path for mod in py_mods}
                assert_equal(files, {os.path.join('subdir', 'pack', '__init__.py'),
                                     os.path.join('subdir', 'a.py'),
                                     os.path.join('subdir', 'b.py')})
            finally:
                os.chdir(current)

    @unittest.skipIf(sys.platform == 'win32', "clean up fails on Windows")
    def test_packages_found(self) -> None:
        current = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            try:
                os.chdir(tmp)
                os.mkdir('pack')
                self.make_file('pack', '__init__.py', content='from . import a, b')
                self.make_file('pack', 'a.py')
                self.make_file('pack', 'b.py')
                opts = parse_options(['-p', 'pack'])
                py_mods, c_mods = collect_build_targets(opts, mypy_options(opts))
                assert_equal(c_mods, [])
                files = {os.path.relpath(mod.path or 'FAIL') for mod in py_mods}
                assert_equal(files, {os.path.join('pack', '__init__.py'),
                                     os.path.join('pack', 'a.py'),
                                     os.path.join('pack', 'b.py')})
            finally:
                os.chdir(current)

    @unittest.skipIf(sys.platform == 'win32', "clean up fails on Windows")
    def test_module_not_found(self) -> None:
        current = os.getcwd()
        captured_output = io.StringIO()
        sys.stdout = captured_output
        with tempfile.TemporaryDirectory() as tmp:
            try:
                os.chdir(tmp)
                self.make_file(tmp, 'mymodule.py', content='import a')
                opts = parse_options(['-m', 'mymodule'])
                py_mods, c_mods = collect_build_targets(opts, mypy_options(opts))
                assert captured_output.getvalue() == ''
            finally:
                sys.stdout = sys.__stdout__
                os.chdir(current)

    def make_file(self, *path: str, content: str = '') -> None:
        file = os.path.join(*path)
        with open(file, 'w') as f:
            f.write(content)

    def run(self, result: Optional[Any] = None) -> Optional[Any]:
        with local_sys_path_set():
            return super().run(result)


class StubgenCliParseSuite(unittest.TestCase):
    def test_walk_packages(self) -> None:
        with ModuleInspect() as m:
            assert_equal(
                set(walk_packages(m, ["mypy.errors"])),
                {"mypy.errors"})

            assert_equal(
                set(walk_packages(m, ["mypy.errors", "mypy.stubgen"])),
                {"mypy.errors", "mypy.stubgen"})

            all_mypy_packages = set(walk_packages(m, ["mypy"]))
            self.assertTrue(all_mypy_packages.issuperset({
                "mypy",
                "mypy.errors",
                "mypy.stubgen",
                "mypy.test",
                "mypy.test.helpers",
            }))


class StubgenUtilSuite(unittest.TestCase):
    """Unit tests for stubgen utility functions."""

    def test_parse_signature(self) -> None:
        self.assert_parse_signature('func()', ('func', [], []))

    def test_parse_signature_with_args(self) -> None:
        self.assert_parse_signature('func(arg)', ('func', ['arg'], []))
        self.assert_parse_signature('do(arg, arg2)', ('do', ['arg', 'arg2'], []))

    def test_parse_signature_with_optional_args(self) -> None:
        self.assert_parse_signature('func([arg])', ('func', [], ['arg']))
        self.assert_parse_signature('func(arg[, arg2])', ('func', ['arg'], ['arg2']))
        self.assert_parse_signature('func([arg[, arg2]])', ('func', [], ['arg', 'arg2']))

    def test_parse_signature_with_default_arg(self) -> None:
        self.assert_parse_signature('func(arg=None)', ('func', [], ['arg']))
        self.assert_parse_signature('func(arg, arg2=None)', ('func', ['arg'], ['arg2']))
        self.assert_parse_signature('func(arg=1, arg2="")', ('func', [], ['arg', 'arg2']))

    def test_parse_signature_with_qualified_function(self) -> None:
        self.assert_parse_signature('ClassName.func(arg)', ('func', ['arg'], []))

    def test_parse_signature_with_kw_only_arg(self) -> None:
        self.assert_parse_signature('ClassName.func(arg, *, arg2=1)',
                                    ('func', ['arg', '*'], ['arg2']))

    def test_parse_signature_with_star_arg(self) -> None:
        self.assert_parse_signature('ClassName.func(arg, *args)',
                                    ('func', ['arg', '*args'], []))

    def test_parse_signature_with_star_star_arg(self) -> None:
        self.assert_parse_signature('ClassName.func(arg, **args)',
                                    ('func', ['arg', '**args'], []))

    def assert_parse_signature(self, sig: str, result: Tuple[str, List[str], List[str]]) -> None:
        assert_equal(parse_signature(sig), result)

    def test_build_signature(self) -> None:
        assert_equal(build_signature([], []), '()')
        assert_equal(build_signature(['arg'], []), '(arg)')
        assert_equal(build_signature(['arg', 'arg2'], []), '(arg, arg2)')
        assert_equal(build_signature(['arg'], ['arg2']), '(arg, arg2=...)')
        assert_equal(build_signature(['arg'], ['arg2', '**x']), '(arg, arg2=..., **x)')

    def test_parse_all_signatures(self) -> None:
        assert_equal(parse_all_signatures(['random text',
                                           '.. function:: fn(arg',
                                           '.. function:: fn()',
                                           '  .. method:: fn2(arg)']),
                     ([('fn', '()'),
                       ('fn2', '(arg)')], []))

    def test_find_unique_signatures(self) -> None:
        assert_equal(find_unique_signatures(
            [('func', '()'),
             ('func', '()'),
             ('func2', '()'),
             ('func2', '(arg)'),
             ('func3', '(arg, arg2)')]),
            [('func', '()'),
             ('func3', '(arg, arg2)')])

    def test_infer_sig_from_docstring(self) -> None:
        assert_equal(infer_sig_from_docstring('\nfunc(x) - y', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x')], ret_type='Any')])
        assert_equal(infer_sig_from_docstring('\nfunc(x)', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x')], ret_type='Any')])

        assert_equal(infer_sig_from_docstring('\nfunc(x, Y_a=None)', 'func'),
                     [FunctionSig(name='func',
                                  args=[ArgSig(name='x'), ArgSig(name='Y_a', default=True)],
                                  ret_type='Any')])

        assert_equal(infer_sig_from_docstring('\nfunc(x, Y_a=3)', 'func'),
                     [FunctionSig(name='func',
                                  args=[ArgSig(name='x'), ArgSig(name='Y_a', default=True)],
                                  ret_type='Any')])

        assert_equal(infer_sig_from_docstring('\nfunc(x, Y_a=[1, 2, 3])', 'func'),
                     [FunctionSig(name='func',
                                  args=[ArgSig(name='x'), ArgSig(name='Y_a', default=True)],
                                  ret_type='Any')])

        assert_equal(infer_sig_from_docstring('\nafunc(x) - y', 'func'), [])
        assert_equal(infer_sig_from_docstring('\nfunc(x, y', 'func'), [])
        assert_equal(infer_sig_from_docstring('\nfunc(x=z(y))', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x', default=True)],
                                  ret_type='Any')])

        assert_equal(infer_sig_from_docstring('\nfunc x', 'func'), [])
        # Try to infer signature from type annotation.
        assert_equal(infer_sig_from_docstring('\nfunc(x: int)', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x', type='int')],
                                  ret_type='Any')])
        assert_equal(infer_sig_from_docstring('\nfunc(x: int=3)', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x', type='int', default=True)],
                                  ret_type='Any')])

        assert_equal(infer_sig_from_docstring('\nfunc(x=3)', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x', type=None, default=True)],
                                  ret_type='Any')])

        assert_equal(infer_sig_from_docstring('\nfunc() -> int', 'func'),
                     [FunctionSig(name='func', args=[], ret_type='int')])

        assert_equal(infer_sig_from_docstring('\nfunc(x: int=3) -> int', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x', type='int', default=True)],
                                  ret_type='int')])

        assert_equal(infer_sig_from_docstring('\nfunc(x: int=3) -> int   \n', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x', type='int', default=True)],
                                  ret_type='int')])

        assert_equal(infer_sig_from_docstring('\nfunc(x: Tuple[int, str]) -> str', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x', type='Tuple[int,str]')],
                                  ret_type='str')])

        assert_equal(
            infer_sig_from_docstring('\nfunc(x: Tuple[int, Tuple[str, int], str], y: int) -> str',
                                     'func'),
            [FunctionSig(name='func',
                         args=[ArgSig(name='x', type='Tuple[int,Tuple[str,int],str]'),
                               ArgSig(name='y', type='int')],
                         ret_type='str')])

        assert_equal(infer_sig_from_docstring('\nfunc(x: foo.bar)', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x', type='foo.bar')],
                                  ret_type='Any')])

        assert_equal(infer_sig_from_docstring('\nfunc(x: list=[1,2,[3,4]])', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x', type='list', default=True)],
                                  ret_type='Any')])

        assert_equal(infer_sig_from_docstring('\nfunc(x: str="nasty[")', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x', type='str', default=True)],
                                  ret_type='Any')])

        assert_equal(infer_sig_from_docstring('\nfunc[(x: foo.bar, invalid]', 'func'), [])

        assert_equal(infer_sig_from_docstring('\nfunc(x: invalid::type<with_template>)', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x', type=None)],
                                  ret_type='Any')])

        assert_equal(infer_sig_from_docstring('\nfunc(x: str="")', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x', type='str', default=True)],
                                  ret_type='Any')])

    def test_infer_sig_from_docstring_duplicate_args(self) -> None:
        assert_equal(infer_sig_from_docstring('\nfunc(x, x) -> str\nfunc(x, y) -> int', 'func'),
                     [FunctionSig(name='func', args=[ArgSig(name='x'), ArgSig(name='y')],
                                  ret_type='int')])

    def test_infer_sig_from_docstring_bad_indentation(self) -> None:
        assert_equal(infer_sig_from_docstring("""
            x
              x
             x
            """, 'func'), None)

    def test_infer_arg_sig_from_anon_docstring(self) -> None:
        assert_equal(infer_arg_sig_from_anon_docstring("(*args, **kwargs)"),
                     [ArgSig(name='*args'), ArgSig(name='**kwargs')])

        assert_equal(
            infer_arg_sig_from_anon_docstring(
                "(x: Tuple[int, Tuple[str, int], str]=(1, ('a', 2), 'y'), y: int=4)"),
            [ArgSig(name='x', type='Tuple[int,Tuple[str,int],str]', default=True),
             ArgSig(name='y', type='int', default=True)])

    def test_infer_prop_type_from_docstring(self) -> None:
        assert_equal(infer_prop_type_from_docstring('str: A string.'), 'str')
        assert_equal(infer_prop_type_from_docstring('Optional[int]: An int.'), 'Optional[int]')
        assert_equal(infer_prop_type_from_docstring('Tuple[int, int]: A tuple.'),
                     'Tuple[int, int]')
        assert_equal(infer_prop_type_from_docstring('\nstr: A string.'), None)

    def test_infer_sig_from_docstring_square_brackets(self) -> None:
        assert infer_sig_from_docstring(
            'fetch_row([maxrows, how]) -- Fetches stuff',
            'fetch_row',
        ) == []

    def test_remove_misplaced_type_comments_1(self) -> None:
        good = """
        \u1234
        def f(x):  # type: (int) -> int

        def g(x):
            # type: (int) -> int

        def h():

            # type: () int

        x = 1  # type: int
        """

        assert_equal(remove_misplaced_type_comments(good), good)

    def test_remove_misplaced_type_comments_2(self) -> None:
        bad = """
        def f(x):
            # type: Callable[[int], int]
            pass

        #  type:  "foo"
        #  type:  'bar'
        x = 1
        # type: int
        """
        bad_fixed = """
        def f(x):

            pass



        x = 1

        """
        assert_equal(remove_misplaced_type_comments(bad), bad_fixed)

    def test_remove_misplaced_type_comments_3(self) -> None:
        bad = '''
        def f(x):
            """docstring"""
            # type: (int) -> int
            pass

        def g(x):
            """docstring
            """
            # type: (int) -> int
            pass
        '''
        bad_fixed = '''
        def f(x):
            """docstring"""

            pass

        def g(x):
            """docstring
            """

            pass
        '''
        assert_equal(remove_misplaced_type_comments(bad), bad_fixed)

    def test_remove_misplaced_type_comments_4(self) -> None:
        bad = """
        def f(x):
            '''docstring'''
            # type: (int) -> int
            pass

        def g(x):
            '''docstring
            '''
            # type: (int) -> int
            pass
        """
        bad_fixed = """
        def f(x):
            '''docstring'''

            pass

        def g(x):
            '''docstring
            '''

            pass
        """
        assert_equal(remove_misplaced_type_comments(bad), bad_fixed)

    def test_remove_misplaced_type_comments_5(self) -> None:
        bad = """
        def f(x):
            # type: (int, List[Any],
            #        float, bool) -> int
            pass

        def g(x):
            # type: (int, List[Any])
            pass
        """
        bad_fixed = """
        def f(x):

            #        float, bool) -> int
            pass

        def g(x):

            pass
        """
        assert_equal(remove_misplaced_type_comments(bad), bad_fixed)

    def test_remove_misplaced_type_comments_bytes(self) -> None:
        original = b"""
        \xbf
        def f(x):  # type: (int) -> int

        def g(x):
            # type: (int) -> int
            pass

        def h():
            # type: int
            pass

        x = 1  # type: int
        """

        dest = b"""
        \xbf
        def f(x):  # type: (int) -> int

        def g(x):
            # type: (int) -> int
            pass

        def h():

            pass

        x = 1  # type: int
        """

        assert_equal(remove_misplaced_type_comments(original), dest)

    @unittest.skipIf(sys.platform == 'win32',
                     'Tests building the paths common ancestor on *nix')
    def test_common_dir_prefix_unix(self) -> None:
        assert common_dir_prefix([]) == '.'
        assert common_dir_prefix(['x.pyi']) == '.'
        assert common_dir_prefix(['./x.pyi']) == '.'
        assert common_dir_prefix(['foo/bar/x.pyi']) == 'foo/bar'
        assert common_dir_prefix(['foo/bar/x.pyi',
                                  'foo/bar/y.pyi']) == 'foo/bar'
        assert common_dir_prefix(['foo/bar/x.pyi', 'foo/y.pyi']) == 'foo'
        assert common_dir_prefix(['foo/x.pyi', 'foo/bar/y.pyi']) == 'foo'
        assert common_dir_prefix(['foo/bar/zar/x.pyi', 'foo/y.pyi']) == 'foo'
        assert common_dir_prefix(['foo/x.pyi', 'foo/bar/zar/y.pyi']) == 'foo'
        assert common_dir_prefix(['foo/bar/zar/x.pyi', 'foo/bar/y.pyi']) == 'foo/bar'
        assert common_dir_prefix(['foo/bar/x.pyi', 'foo/bar/zar/y.pyi']) == 'foo/bar'
        assert common_dir_prefix([r'foo/bar\x.pyi']) == 'foo'
        assert common_dir_prefix([r'foo\bar/x.pyi']) == r'foo\bar'

    @unittest.skipIf(sys.platform != 'win32',
                     'Tests building the paths common ancestor on Windows')
    def test_common_dir_prefix_win(self) -> None:
        assert common_dir_prefix(['x.pyi']) == '.'
        assert common_dir_prefix([r'.\x.pyi']) == '.'
        assert common_dir_prefix([r'foo\bar\x.pyi']) == r'foo\bar'
        assert common_dir_prefix([r'foo\bar\x.pyi',
                                  r'foo\bar\y.pyi']) == r'foo\bar'
        assert common_dir_prefix([r'foo\bar\x.pyi', r'foo\y.pyi']) == 'foo'
        assert common_dir_prefix([r'foo\x.pyi', r'foo\bar\y.pyi']) == 'foo'
        assert common_dir_prefix([r'foo\bar\zar\x.pyi', r'foo\y.pyi']) == 'foo'
        assert common_dir_prefix([r'foo\x.pyi', r'foo\bar\zar\y.pyi']) == 'foo'
        assert common_dir_prefix([r'foo\bar\zar\x.pyi', r'foo\bar\y.pyi']) == r'foo\bar'
        assert common_dir_prefix([r'foo\bar\x.pyi', r'foo\bar\zar\y.pyi']) == r'foo\bar'
        assert common_dir_prefix([r'foo/bar\x.pyi']) == r'foo\bar'
        assert common_dir_prefix([r'foo\bar/x.pyi']) == r'foo\bar'
        assert common_dir_prefix([r'foo/bar/x.pyi']) == r'foo\bar'


class StubgenHelpersSuite(unittest.TestCase):
    def test_is_blacklisted_path(self) -> None:
        assert not is_blacklisted_path('foo/bar.py')
        assert not is_blacklisted_path('foo.py')
        assert not is_blacklisted_path('foo/xvendor/bar.py')
        assert not is_blacklisted_path('foo/vendorx/bar.py')
        assert is_blacklisted_path('foo/vendor/bar.py')
        assert is_blacklisted_path('foo/vendored/bar.py')
        assert is_blacklisted_path('foo/vendored/bar/thing.py')
        assert is_blacklisted_path('foo/six.py')

    def test_is_non_library_module(self) -> None:
        assert not is_non_library_module('foo')
        assert not is_non_library_module('foo.bar')

        # The following could be test modules, but we are very conservative and
        # don't treat them as such since they could plausibly be real modules.
        assert not is_non_library_module('foo.bartest')
        assert not is_non_library_module('foo.bartests')
        assert not is_non_library_module('foo.testbar')

        assert is_non_library_module('foo.test')
        assert is_non_library_module('foo.test.foo')
        assert is_non_library_module('foo.tests')
        assert is_non_library_module('foo.tests.foo')
        assert is_non_library_module('foo.testing.foo')
        assert is_non_library_module('foo.SelfTest.foo')

        assert is_non_library_module('foo.test_bar')
        assert is_non_library_module('foo.bar_tests')
        assert is_non_library_module('foo.testing')
        assert is_non_library_module('foo.conftest')
        assert is_non_library_module('foo.bar_test_util')
        assert is_non_library_module('foo.bar_test_utils')
        assert is_non_library_module('foo.bar_test_base')

        assert is_non_library_module('foo.setup')

        assert is_non_library_module('foo.__main__')


class StubgenPythonSuite(DataSuite):
    """Data-driven end-to-end test cases that generate stub files.

    You can use these magic test case name suffixes:

    *_semanal
        Run semantic analysis (slow as this uses real stubs -- only use
        when necessary)
    *_import
        Import module and perform runtime introspection (in the current
        process!)

    You can use these magic comments:

    # flags: --some-stubgen-option ...
        Specify custom stubgen options

    # modules: module1 module2 ...
        Specify which modules to output (by default only 'main')
    """

    required_out_section = True
    base_path = '.'
    files = ['stubgen.test']

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        with local_sys_path_set():
            self.run_case_inner(testcase)

    def run_case_inner(self, testcase: DataDrivenTestCase) -> None:
        extra = []  # Extra command-line args
        mods = []  # Module names to process
        source = '\n'.join(testcase.input)
        for file, content in testcase.files + [('./main.py', source)]:
            # Strip ./ prefix and .py suffix.
            mod = file[2:-3].replace('/', '.')
            if mod.endswith('.__init__'):
                mod, _, _ = mod.rpartition('.')
            mods.append(mod)
            if '-p ' not in source:
                extra.extend(['-m', mod])
            with open(file, 'w') as f:
                f.write(content)

        options = self.parse_flags(source, extra)
        modules = self.parse_modules(source)
        out_dir = 'out'
        try:
            try:
                if not testcase.name.endswith('_import'):
                    options.no_import = True
                if not testcase.name.endswith('_semanal'):
                    options.parse_only = True
                generate_stubs(options)
                a: List[str] = []
                for module in modules:
                    fnam = module_to_path(out_dir, module)
                    self.add_file(fnam, a, header=len(modules) > 1)
            except CompileError as e:
                a = e.messages
            assert_string_arrays_equal(testcase.output, a,
                                       'Invalid output ({}, line {})'.format(
                                           testcase.file, testcase.line))
        finally:
            for mod in mods:
                if mod in sys.modules:
                    del sys.modules[mod]
            shutil.rmtree(out_dir)

    def parse_flags(self, program_text: str, extra: List[str]) -> Options:
        flags = re.search('# flags: (.*)$', program_text, flags=re.MULTILINE)
        if flags:
            flag_list = flags.group(1).split()
        else:
            flag_list = []
        options = parse_options(flag_list + extra)
        if '--verbose' not in flag_list:
            options.quiet = True
        else:
            options.verbose = True
        return options

    def parse_modules(self, program_text: str) -> List[str]:
        modules = re.search('# modules: (.*)$', program_text, flags=re.MULTILINE)
        if modules:
            return modules.group(1).split()
        else:
            return ['main']

    def add_file(self, path: str, result: List[str], header: bool) -> None:
        if not os.path.exists(path):
            result.append('<%s was not generated>' % path.replace('\\', '/'))
            return
        if header:
            result.append('# {}'.format(path[4:]))
        with open(path, encoding='utf8') as file:
            result.extend(file.read().splitlines())


self_arg = ArgSig(name='self')


class TestBaseClass:
    pass


class TestClass(TestBaseClass):
    pass


class StubgencSuite(unittest.TestCase):
    """Unit tests for stub generation from C modules using introspection.

    Note that these don't cover a lot!
    """

    def test_infer_hash_sig(self) -> None:
        assert_equal(infer_method_sig('__hash__'), [self_arg])

    def test_infer_getitem_sig(self) -> None:
        assert_equal(infer_method_sig('__getitem__'), [self_arg, ArgSig(name='index')])

    def test_infer_setitem_sig(self) -> None:
        assert_equal(infer_method_sig('__setitem__'),
                     [self_arg, ArgSig(name='index'), ArgSig(name='object')])

    def test_infer_binary_op_sig(self) -> None:
        for op in ('eq', 'ne', 'lt', 'le', 'gt', 'ge',
                   'add', 'radd', 'sub', 'rsub', 'mul', 'rmul'):
            assert_equal(infer_method_sig('__%s__' % op), [self_arg, ArgSig(name='other')])

    def test_infer_unary_op_sig(self) -> None:
        for op in ('neg', 'pos'):
            assert_equal(infer_method_sig('__%s__' % op), [self_arg])

    def test_generate_c_type_stub_no_crash_for_object(self) -> None:
        output: List[str] = []
        mod = ModuleType('module', '')  # any module is fine
        imports: List[str] = []
        generate_c_type_stub(mod, 'alias', object, output, imports)
        assert_equal(imports, [])
        assert_equal(output[0], 'class alias:')

    def test_generate_c_type_stub_variable_type_annotation(self) -> None:
        # This class mimics the stubgen unit test 'testClassVariable'
        class TestClassVariableCls:
            x = 1

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType('module', '')  # any module is fine
        generate_c_type_stub(mod, 'C', TestClassVariableCls, output, imports)
        assert_equal(imports, [])
        assert_equal(output, ['class C:', '    x: ClassVar[int] = ...'])

    def test_generate_c_type_inheritance(self) -> None:
        class TestClass(KeyError):
            pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType('module, ')
        generate_c_type_stub(mod, 'C', TestClass, output, imports)
        assert_equal(output, ['class C(KeyError): ...', ])
        assert_equal(imports, [])

    def test_generate_c_type_inheritance_same_module(self) -> None:
        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType(TestBaseClass.__module__, '')
        generate_c_type_stub(mod, 'C', TestClass, output, imports)
        assert_equal(output, ['class C(TestBaseClass): ...', ])
        assert_equal(imports, [])

    def test_generate_c_type_inheritance_other_module(self) -> None:
        import argparse

        class TestClass(argparse.Action):
            pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType('module', '')
        generate_c_type_stub(mod, 'C', TestClass, output, imports)
        assert_equal(output, ['class C(argparse.Action): ...', ])
        assert_equal(imports, ['import argparse'])

    def test_generate_c_type_inheritance_builtin_type(self) -> None:
        class TestClass(type):
            pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType('module', '')
        generate_c_type_stub(mod, 'C', TestClass, output, imports)
        assert_equal(output, ['class C(type): ...', ])
        assert_equal(imports, [])

    def test_generate_c_type_with_docstring(self) -> None:
        class TestClass:
            def test(self, arg0: str) -> None:
                """
                test(self: TestClass, arg0: int)
                """
                pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType(TestClass.__module__, '')
        generate_c_function_stub(mod, 'test', TestClass.test, output, imports,
                                 self_var='self', class_name='TestClass')
        assert_equal(output, ['def test(self, arg0: int) -> Any: ...'])
        assert_equal(imports, [])

    def test_generate_c_type_with_docstring_no_self_arg(self) -> None:
        class TestClass:
            def test(self, arg0: str) -> None:
                """
                test(arg0: int)
                """
                pass
        output = []  # type: List[str]
        imports = []  # type: List[str]
        mod = ModuleType(TestClass.__module__, '')
        generate_c_function_stub(mod, 'test', TestClass.test, output, imports,
                                 self_var='self', class_name='TestClass')
        assert_equal(output, ['def test(self, arg0: int) -> Any: ...'])
        assert_equal(imports, [])

    def test_generate_c_type_classmethod(self) -> None:
        class TestClass:
            @classmethod
            def test(cls, arg0: str) -> None:
                pass
        output = []  # type: List[str]
        imports = []  # type: List[str]
        mod = ModuleType(TestClass.__module__, '')
        generate_c_function_stub(mod, 'test', TestClass.test, output, imports,
                                 self_var='cls', class_name='TestClass')
        assert_equal(output, ['def test(cls, *args, **kwargs) -> Any: ...'])
        assert_equal(imports, [])

    def test_generate_c_type_with_docstring_empty_default(self) -> None:
        class TestClass:
            def test(self, arg0: str = "") -> None:
                """
                test(self: TestClass, arg0: str = "")
                """
                pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType(TestClass.__module__, '')
        generate_c_function_stub(mod, 'test', TestClass.test, output, imports,
                                 self_var='self', class_name='TestClass')
        assert_equal(output, ['def test(self, arg0: str = ...) -> Any: ...'])
        assert_equal(imports, [])

    def test_generate_c_function_other_module_arg(self) -> None:
        """Test that if argument references type from other module, module will be imported."""
        # Provide different type in python spec than in docstring to make sure, that docstring
        # information is used.
        def test(arg0: str) -> None:
            """
            test(arg0: argparse.Action)
            """
            pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType(self.__module__, '')
        generate_c_function_stub(mod, 'test', test, output, imports)
        assert_equal(output, ['def test(arg0: argparse.Action) -> Any: ...'])
        assert_equal(imports, ['import argparse'])

    def test_generate_c_function_same_module_arg(self) -> None:
        """Test that if argument references type from same module but using full path, no module
        will be imported, and type specification will be striped to local reference.
        """
        # Provide different type in python spec than in docstring to make sure, that docstring
        # information is used.
        def test(arg0: str) -> None:
            """
            test(arg0: argparse.Action)
            """
            pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType('argparse', '')
        generate_c_function_stub(mod, 'test', test, output, imports)
        assert_equal(output, ['def test(arg0: Action) -> Any: ...'])
        assert_equal(imports, [])

    def test_generate_c_function_other_module_ret(self) -> None:
        """Test that if return type references type from other module, module will be imported."""
        def test(arg0: str) -> None:
            """
            test(arg0: str) -> argparse.Action
            """
            pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType(self.__module__, '')
        generate_c_function_stub(mod, 'test', test, output, imports)
        assert_equal(output, ['def test(arg0: str) -> argparse.Action: ...'])
        assert_equal(imports, ['import argparse'])

    def test_generate_c_function_same_module_ret(self) -> None:
        """Test that if return type references type from same module but using full path,
        no module will be imported, and type specification will be striped to local reference.
        """
        def test(arg0: str) -> None:
            """
            test(arg0: str) -> argparse.Action
            """
            pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType('argparse', '')
        generate_c_function_stub(mod, 'test', test, output, imports)
        assert_equal(output, ['def test(arg0: str) -> Action: ...'])
        assert_equal(imports, [])

    def test_generate_c_property_with_pybind11(self) -> None:
        """Signatures included by PyBind11 inside property.fget are read."""
        class TestClass:
            def get_attribute(self) -> None:
                """
                (self: TestClass) -> str
                """
                pass
            attribute = property(get_attribute, doc="")

        output: List[str] = []
        generate_c_property_stub('attribute', TestClass.attribute, [], [], output, readonly=True)
        assert_equal(output, ['@property', 'def attribute(self) -> str: ...'])

    def test_generate_c_type_with_single_arg_generic(self) -> None:
        class TestClass:
            def test(self, arg0: str) -> None:
                """
                test(self: TestClass, arg0: List[int])
                """
                pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType(TestClass.__module__, '')
        generate_c_function_stub(mod, 'test', TestClass.test, output, imports,
                                 self_var='self', class_name='TestClass')
        assert_equal(output, ['def test(self, arg0: List[int]) -> Any: ...'])
        assert_equal(imports, [])

    def test_generate_c_type_with_double_arg_generic(self) -> None:
        class TestClass:
            def test(self, arg0: str) -> None:
                """
                test(self: TestClass, arg0: Dict[str, int])
                """
                pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType(TestClass.__module__, '')
        generate_c_function_stub(mod, 'test', TestClass.test, output, imports,
                                 self_var='self', class_name='TestClass')
        assert_equal(output, ['def test(self, arg0: Dict[str,int]) -> Any: ...'])
        assert_equal(imports, [])

    def test_generate_c_type_with_nested_generic(self) -> None:
        class TestClass:
            def test(self, arg0: str) -> None:
                """
                test(self: TestClass, arg0: Dict[str, List[int]])
                """
                pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType(TestClass.__module__, '')
        generate_c_function_stub(mod, 'test', TestClass.test, output, imports,
                                 self_var='self', class_name='TestClass')
        assert_equal(output, ['def test(self, arg0: Dict[str,List[int]]) -> Any: ...'])
        assert_equal(imports, [])

    def test_generate_c_type_with_generic_using_other_module_first(self) -> None:
        class TestClass:
            def test(self, arg0: str) -> None:
                """
                test(self: TestClass, arg0: Dict[argparse.Action, int])
                """
                pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType(TestClass.__module__, '')
        generate_c_function_stub(mod, 'test', TestClass.test, output, imports,
                                 self_var='self', class_name='TestClass')
        assert_equal(output, ['def test(self, arg0: Dict[argparse.Action,int]) -> Any: ...'])
        assert_equal(imports, ['import argparse'])

    def test_generate_c_type_with_generic_using_other_module_last(self) -> None:
        class TestClass:
            def test(self, arg0: str) -> None:
                """
                test(self: TestClass, arg0: Dict[str, argparse.Action])
                """
                pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType(TestClass.__module__, '')
        generate_c_function_stub(mod, 'test', TestClass.test, output, imports,
                                 self_var='self', class_name='TestClass')
        assert_equal(output, ['def test(self, arg0: Dict[str,argparse.Action]) -> Any: ...'])
        assert_equal(imports, ['import argparse'])

    def test_generate_c_type_with_overload_pybind11(self) -> None:
        class TestClass:
            def __init__(self, arg0: str) -> None:
                """
                __init__(*args, **kwargs)
                Overloaded function.

                1. __init__(self: TestClass, arg0: str) -> None

                2. __init__(self: TestClass, arg0: str, arg1: str) -> None
                """
                pass

        output: List[str] = []
        imports: List[str] = []
        mod = ModuleType(TestClass.__module__, '')
        generate_c_function_stub(mod, '__init__', TestClass.__init__, output, imports,
                                 self_var='self', class_name='TestClass')
        assert_equal(output, [
            '@overload',
            'def __init__(self, arg0: str) -> None: ...',
            '@overload',
            'def __init__(self, arg0: str, arg1: str) -> None: ...',
            '@overload',
            'def __init__(*args, **kwargs) -> Any: ...'])
        assert_equal(set(imports), {'from typing import overload'})


class ArgSigSuite(unittest.TestCase):
    def test_repr(self) -> None:
        assert_equal(repr(ArgSig(name='asd"dsa')),
                     "ArgSig(name='asd\"dsa', type=None, default=False)")
        assert_equal(repr(ArgSig(name="asd'dsa")),
                     'ArgSig(name="asd\'dsa", type=None, default=False)')
        assert_equal(repr(ArgSig("func", 'str')),
                     "ArgSig(name='func', type='str', default=False)")
        assert_equal(repr(ArgSig("func", 'str', default=True)),
                     "ArgSig(name='func', type='str', default=True)")


class IsValidTypeSuite(unittest.TestCase):
    def test_is_valid_type(self) -> None:
        assert is_valid_type('int')
        assert is_valid_type('str')
        assert is_valid_type('Foo_Bar234')
        assert is_valid_type('foo.bar')
        assert is_valid_type('List[int]')
        assert is_valid_type('Dict[str, int]')
        assert is_valid_type('None')
        assert not is_valid_type('foo-bar')
        assert not is_valid_type('x->y')
        assert not is_valid_type('True')
        assert not is_valid_type('False')
        assert not is_valid_type('x,y')
        assert not is_valid_type('x, y')


class ModuleInspectSuite(unittest.TestCase):
    def test_python_module(self) -> None:
        with ModuleInspect() as m:
            p = m.get_package_properties('inspect')
            assert p is not None
            assert p.name == 'inspect'
            assert p.file
            assert p.path is None
            assert p.is_c_module is False
            assert p.subpackages == []

    def test_python_package(self) -> None:
        with ModuleInspect() as m:
            p = m.get_package_properties('unittest')
            assert p is not None
            assert p.name == 'unittest'
            assert p.file
            assert p.path
            assert p.is_c_module is False
            assert p.subpackages
            assert all(sub.startswith('unittest.') for sub in p.subpackages)

    def test_c_module(self) -> None:
        with ModuleInspect() as m:
            p = m.get_package_properties('_socket')
            assert p is not None
            assert p.name == '_socket'
            assert p.path is None
            assert p.is_c_module is True
            assert p.subpackages == []

    def test_non_existent(self) -> None:
        with ModuleInspect() as m:
            with self.assertRaises(InspectError) as e:
                m.get_package_properties('foobar-non-existent')
            assert str(e.exception) == "No module named 'foobar-non-existent'"


def module_to_path(out_dir: str, module: str) -> str:
    fnam = os.path.join(out_dir, '{}.pyi'.format(module.replace('.', '/')))
    if not os.path.exists(fnam):
        alt_fnam = fnam.replace('.pyi', '/__init__.pyi')
        if os.path.exists(alt_fnam):
            return alt_fnam
    return fnam
