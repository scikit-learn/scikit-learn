# -*- coding: utf-8 -*-
"""Tests for the key interactiveshell module.

Historically the main classes in interactiveshell have been under-tested.  This
module should grow as many single-method tests as possible to trap many of the
recurring bugs we seem to encounter with high-level interaction.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import ast
import os
import signal
import shutil
import sys
import tempfile
import unittest
from unittest import mock

from os.path import join

import nose.tools as nt

from IPython.core.error import InputRejected
from IPython.core.inputtransformer import InputTransformer
from IPython.testing.decorators import (
    skipif, skip_win32, onlyif_unicode_paths, onlyif_cmds_exist,
)
from IPython.testing import tools as tt
from IPython.utils.process import find_cmd

#-----------------------------------------------------------------------------
# Globals
#-----------------------------------------------------------------------------
# This is used by every single test, no point repeating it ad nauseam
ip = get_ipython()

#-----------------------------------------------------------------------------
# Tests
#-----------------------------------------------------------------------------

class DerivedInterrupt(KeyboardInterrupt):
    pass

class InteractiveShellTestCase(unittest.TestCase):
    def test_naked_string_cells(self):
        """Test that cells with only naked strings are fully executed"""
        # First, single-line inputs
        ip.run_cell('"a"\n')
        self.assertEqual(ip.user_ns['_'], 'a')
        # And also multi-line cells
        ip.run_cell('"""a\nb"""\n')
        self.assertEqual(ip.user_ns['_'], 'a\nb')

    def test_run_empty_cell(self):
        """Just make sure we don't get a horrible error with a blank
        cell of input. Yes, I did overlook that."""
        old_xc = ip.execution_count
        res = ip.run_cell('')
        self.assertEqual(ip.execution_count, old_xc)
        self.assertEqual(res.execution_count, None)

    def test_run_cell_multiline(self):
        """Multi-block, multi-line cells must execute correctly.
        """
        src = '\n'.join(["x=1",
                         "y=2",
                         "if 1:",
                         "    x += 1",
                         "    y += 1",])
        res = ip.run_cell(src)
        self.assertEqual(ip.user_ns['x'], 2)
        self.assertEqual(ip.user_ns['y'], 3)
        self.assertEqual(res.success, True)
        self.assertEqual(res.result, None)

    def test_multiline_string_cells(self):
        "Code sprinkled with multiline strings should execute (GH-306)"
        ip.run_cell('tmp=0')
        self.assertEqual(ip.user_ns['tmp'], 0)
        res = ip.run_cell('tmp=1;"""a\nb"""\n')
        self.assertEqual(ip.user_ns['tmp'], 1)
        self.assertEqual(res.success, True)
        self.assertEqual(res.result, "a\nb")

    def test_dont_cache_with_semicolon(self):
        "Ending a line with semicolon should not cache the returned object (GH-307)"
        oldlen = len(ip.user_ns['Out'])
        for cell in ['1;', '1;1;']:
            res = ip.run_cell(cell, store_history=True)
            newlen = len(ip.user_ns['Out'])
            self.assertEqual(oldlen, newlen)
            self.assertIsNone(res.result)
        i = 0
        #also test the default caching behavior
        for cell in ['1', '1;1']:
            ip.run_cell(cell, store_history=True)
            newlen = len(ip.user_ns['Out'])
            i += 1
            self.assertEqual(oldlen+i, newlen)

    def test_syntax_error(self):
        res = ip.run_cell("raise = 3")
        self.assertIsInstance(res.error_before_exec, SyntaxError)

    def test_In_variable(self):
        "Verify that In variable grows with user input (GH-284)"
        oldlen = len(ip.user_ns['In'])
        ip.run_cell('1;', store_history=True)
        newlen = len(ip.user_ns['In'])
        self.assertEqual(oldlen+1, newlen)
        self.assertEqual(ip.user_ns['In'][-1],'1;')
        
    def test_magic_names_in_string(self):
        ip.run_cell('a = """\n%exit\n"""')
        self.assertEqual(ip.user_ns['a'], '\n%exit\n')
    
    def test_trailing_newline(self):
        """test that running !(command) does not raise a SyntaxError"""
        ip.run_cell('!(true)\n', False)
        ip.run_cell('!(true)\n\n\n', False)
    
    def test_gh_597(self):
        """Pretty-printing lists of objects with non-ascii reprs may cause
        problems."""
        class Spam(object):
          def __repr__(self):
            return "\xe9"*50
        import IPython.core.formatters
        f = IPython.core.formatters.PlainTextFormatter()
        f([Spam(),Spam()])
    

    def test_future_flags(self):
        """Check that future flags are used for parsing code (gh-777)"""
        ip.run_cell('from __future__ import barry_as_FLUFL')
        try:
            ip.run_cell('prfunc_return_val = 1 <> 2')
            assert 'prfunc_return_val' in ip.user_ns
        finally:
            # Reset compiler flags so we don't mess up other tests.
            ip.compile.reset_compiler_flags()

    def test_can_pickle(self):
        "Can we pickle objects defined interactively (GH-29)"
        ip = get_ipython()
        ip.reset()
        ip.run_cell(("class Mylist(list):\n"
                     "    def __init__(self,x=[]):\n"
                     "        list.__init__(self,x)"))
        ip.run_cell("w=Mylist([1,2,3])")
        
        from pickle import dumps
        
        # We need to swap in our main module - this is only necessary
        # inside the test framework, because IPython puts the interactive module
        # in place (but the test framework undoes this).
        _main = sys.modules['__main__']
        sys.modules['__main__'] = ip.user_module
        try:
            res = dumps(ip.user_ns["w"])
        finally:
            sys.modules['__main__'] = _main
        self.assertTrue(isinstance(res, bytes))
        
    def test_global_ns(self):
        "Code in functions must be able to access variables outside them."
        ip = get_ipython()
        ip.run_cell("a = 10")
        ip.run_cell(("def f(x):\n"
                     "    return x + a"))
        ip.run_cell("b = f(12)")
        self.assertEqual(ip.user_ns["b"], 22)

    def test_bad_custom_tb(self):
        """Check that InteractiveShell is protected from bad custom exception handlers"""
        ip.set_custom_exc((IOError,), lambda etype,value,tb: 1/0)
        self.assertEqual(ip.custom_exceptions, (IOError,))
        with tt.AssertPrints("Custom TB Handler failed", channel='stderr'):
            ip.run_cell(u'raise IOError("foo")')
        self.assertEqual(ip.custom_exceptions, ())

    def test_bad_custom_tb_return(self):
        """Check that InteractiveShell is protected from bad return types in custom exception handlers"""
        ip.set_custom_exc((NameError,),lambda etype,value,tb, tb_offset=None: 1)
        self.assertEqual(ip.custom_exceptions, (NameError,))
        with tt.AssertPrints("Custom TB Handler failed", channel='stderr'):
            ip.run_cell(u'a=abracadabra')
        self.assertEqual(ip.custom_exceptions, ())

    def test_drop_by_id(self):
        myvars = {"a":object(), "b":object(), "c": object()}
        ip.push(myvars, interactive=False)
        for name in myvars:
            assert name in ip.user_ns, name
            assert name in ip.user_ns_hidden, name
        ip.user_ns['b'] = 12
        ip.drop_by_id(myvars)
        for name in ["a", "c"]:
            assert name not in ip.user_ns, name
            assert name not in ip.user_ns_hidden, name
        assert ip.user_ns['b'] == 12
        ip.reset()

    def test_var_expand(self):
        ip.user_ns['f'] = u'Ca\xf1o'
        self.assertEqual(ip.var_expand(u'echo $f'), u'echo Ca\xf1o')
        self.assertEqual(ip.var_expand(u'echo {f}'), u'echo Ca\xf1o')
        self.assertEqual(ip.var_expand(u'echo {f[:-1]}'), u'echo Ca\xf1')
        self.assertEqual(ip.var_expand(u'echo {1*2}'), u'echo 2')
        
        self.assertEqual(ip.var_expand(u"grep x | awk '{print $1}'"), u"grep x | awk '{print $1}'")

        ip.user_ns['f'] = b'Ca\xc3\xb1o'
        # This should not raise any exception:
        ip.var_expand(u'echo $f')
   
    def test_var_expand_local(self):
        """Test local variable expansion in !system and %magic calls"""
        # !system
        ip.run_cell('def test():\n'
                    '    lvar = "ttt"\n'
                    '    ret = !echo {lvar}\n'
                    '    return ret[0]\n')
        res = ip.user_ns['test']()
        nt.assert_in('ttt', res)
        
        # %magic
        ip.run_cell('def makemacro():\n'
                    '    macroname = "macro_var_expand_locals"\n'
                    '    %macro {macroname} codestr\n')
        ip.user_ns['codestr'] = "str(12)"
        ip.run_cell('makemacro()')
        nt.assert_in('macro_var_expand_locals', ip.user_ns)
    
    def test_var_expand_self(self):
        """Test variable expansion with the name 'self', which was failing.
        
        See https://github.com/ipython/ipython/issues/1878#issuecomment-7698218
        """
        ip.run_cell('class cTest:\n'
                    '  classvar="see me"\n'
                    '  def test(self):\n'
                    '    res = !echo Variable: {self.classvar}\n'
                    '    return res[0]\n')
        nt.assert_in('see me', ip.user_ns['cTest']().test())
    
    def test_bad_var_expand(self):
        """var_expand on invalid formats shouldn't raise"""
        # SyntaxError
        self.assertEqual(ip.var_expand(u"{'a':5}"), u"{'a':5}")
        # NameError
        self.assertEqual(ip.var_expand(u"{asdf}"), u"{asdf}")
        # ZeroDivisionError
        self.assertEqual(ip.var_expand(u"{1/0}"), u"{1/0}")
    
    def test_silent_postexec(self):
        """run_cell(silent=True) doesn't invoke pre/post_run_cell callbacks"""
        pre_explicit = mock.Mock()
        pre_always = mock.Mock()
        post_explicit = mock.Mock()
        post_always = mock.Mock()
        all_mocks = [pre_explicit, pre_always, post_explicit, post_always]
        
        ip.events.register('pre_run_cell', pre_explicit)
        ip.events.register('pre_execute', pre_always)
        ip.events.register('post_run_cell', post_explicit)
        ip.events.register('post_execute', post_always)
        
        try:
            ip.run_cell("1", silent=True)
            assert pre_always.called
            assert not pre_explicit.called
            assert post_always.called
            assert not post_explicit.called
            # double-check that non-silent exec did what we expected
            # silent to avoid
            ip.run_cell("1")
            assert pre_explicit.called
            assert post_explicit.called
            info, = pre_explicit.call_args[0]
            result, = post_explicit.call_args[0]
            self.assertEqual(info, result.info)
            # check that post hooks are always called
            [m.reset_mock() for m in all_mocks]
            ip.run_cell("syntax error")
            assert pre_always.called
            assert pre_explicit.called
            assert post_always.called
            assert post_explicit.called
            info, = pre_explicit.call_args[0]
            result, = post_explicit.call_args[0]
            self.assertEqual(info, result.info)
        finally:
            # remove post-exec
            ip.events.unregister('pre_run_cell', pre_explicit)
            ip.events.unregister('pre_execute', pre_always)
            ip.events.unregister('post_run_cell', post_explicit)
            ip.events.unregister('post_execute', post_always)
    
    def test_silent_noadvance(self):
        """run_cell(silent=True) doesn't advance execution_count"""
        ec = ip.execution_count
        # silent should force store_history=False
        ip.run_cell("1", store_history=True, silent=True)
        
        self.assertEqual(ec, ip.execution_count)
        # double-check that non-silent exec did what we expected
        # silent to avoid
        ip.run_cell("1", store_history=True)
        self.assertEqual(ec+1, ip.execution_count)
    
    def test_silent_nodisplayhook(self):
        """run_cell(silent=True) doesn't trigger displayhook"""
        d = dict(called=False)
        
        trap = ip.display_trap
        save_hook = trap.hook
        
        def failing_hook(*args, **kwargs):
            d['called'] = True
        
        try:
            trap.hook = failing_hook
            res = ip.run_cell("1", silent=True)
            self.assertFalse(d['called'])
            self.assertIsNone(res.result)
            # double-check that non-silent exec did what we expected
            # silent to avoid
            ip.run_cell("1")
            self.assertTrue(d['called'])
        finally:
            trap.hook = save_hook

    def test_ofind_line_magic(self):
        from IPython.core.magic import register_line_magic
        
        @register_line_magic
        def lmagic(line):
            "A line magic"

        # Get info on line magic
        lfind = ip._ofind('lmagic')
        info = dict(found=True, isalias=False, ismagic=True,
                    namespace = 'IPython internal', obj= lmagic.__wrapped__,
                    parent = None)
        nt.assert_equal(lfind, info)
        
    def test_ofind_cell_magic(self):
        from IPython.core.magic import register_cell_magic
        
        @register_cell_magic
        def cmagic(line, cell):
            "A cell magic"

        # Get info on cell magic
        find = ip._ofind('cmagic')
        info = dict(found=True, isalias=False, ismagic=True,
                    namespace = 'IPython internal', obj= cmagic.__wrapped__,
                    parent = None)
        nt.assert_equal(find, info)

    def test_ofind_property_with_error(self):
        class A(object):
            @property
            def foo(self):
                raise NotImplementedError()
        a = A()

        found = ip._ofind('a.foo', [('locals', locals())])
        info = dict(found=True, isalias=False, ismagic=False,
                    namespace='locals', obj=A.foo, parent=a)
        nt.assert_equal(found, info)

    def test_ofind_multiple_attribute_lookups(self):
        class A(object):
            @property
            def foo(self):
                raise NotImplementedError()

        a = A()
        a.a = A()
        a.a.a = A()

        found = ip._ofind('a.a.a.foo', [('locals', locals())])
        info = dict(found=True, isalias=False, ismagic=False,
                    namespace='locals', obj=A.foo, parent=a.a.a)
        nt.assert_equal(found, info)

    def test_ofind_slotted_attributes(self):
        class A(object):
            __slots__ = ['foo']
            def __init__(self):
                self.foo = 'bar'

        a = A()
        found = ip._ofind('a.foo', [('locals', locals())])
        info = dict(found=True, isalias=False, ismagic=False,
                    namespace='locals', obj=a.foo, parent=a)
        nt.assert_equal(found, info)

        found = ip._ofind('a.bar', [('locals', locals())])
        info = dict(found=False, isalias=False, ismagic=False,
                    namespace=None, obj=None, parent=a)
        nt.assert_equal(found, info)

    def test_ofind_prefers_property_to_instance_level_attribute(self):
        class A(object):
            @property
            def foo(self):
                return 'bar'
        a = A()
        a.__dict__['foo'] = 'baz'
        nt.assert_equal(a.foo, 'bar')
        found = ip._ofind('a.foo', [('locals', locals())])
        nt.assert_is(found['obj'], A.foo)

    def test_custom_syntaxerror_exception(self):
        called = []
        def my_handler(shell, etype, value, tb, tb_offset=None):
            called.append(etype)
            shell.showtraceback((etype, value, tb), tb_offset=tb_offset)

        ip.set_custom_exc((SyntaxError,), my_handler)
        try:
            ip.run_cell("1f")
            # Check that this was called, and only once.
            self.assertEqual(called, [SyntaxError])
        finally:
            # Reset the custom exception hook
            ip.set_custom_exc((), None)

    def test_custom_exception(self):
        called = []
        def my_handler(shell, etype, value, tb, tb_offset=None):
            called.append(etype)
            shell.showtraceback((etype, value, tb), tb_offset=tb_offset)
        
        ip.set_custom_exc((ValueError,), my_handler)
        try:
            res = ip.run_cell("raise ValueError('test')")
            # Check that this was called, and only once.
            self.assertEqual(called, [ValueError])
            # Check that the error is on the result object
            self.assertIsInstance(res.error_in_exec, ValueError)
        finally:
            # Reset the custom exception hook
            ip.set_custom_exc((), None)
    
    def test_mktempfile(self):
        filename = ip.mktempfile()
        # Check that we can open the file again on Windows
        with open(filename, 'w') as f:
            f.write('abc')

        filename = ip.mktempfile(data='blah')
        with open(filename, 'r') as f:
            self.assertEqual(f.read(), 'blah')

    def test_new_main_mod(self):
        # Smoketest to check that this accepts a unicode module name
        name = u'jiefmw'
        mod = ip.new_main_mod(u'%s.py' % name, name)
        self.assertEqual(mod.__name__, name)

    def test_get_exception_only(self):
        try:
            raise KeyboardInterrupt
        except KeyboardInterrupt:
            msg = ip.get_exception_only()
        self.assertEqual(msg, 'KeyboardInterrupt\n')

        try:
            raise DerivedInterrupt("foo")
        except KeyboardInterrupt:
            msg = ip.get_exception_only()
        self.assertEqual(msg, 'IPython.core.tests.test_interactiveshell.DerivedInterrupt: foo\n')

    def test_inspect_text(self):
        ip.run_cell('a = 5')
        text = ip.object_inspect_text('a')
        self.assertIsInstance(text, str)

    def test_last_execution_result(self):
        """ Check that last execution result gets set correctly (GH-10702) """
        result = ip.run_cell('a = 5; a')
        self.assertTrue(ip.last_execution_succeeded)
        self.assertEqual(ip.last_execution_result.result, 5)

        result = ip.run_cell('a = x_invalid_id_x')
        self.assertFalse(ip.last_execution_succeeded)
        self.assertFalse(ip.last_execution_result.success)
        self.assertIsInstance(ip.last_execution_result.error_in_exec, NameError)


class TestSafeExecfileNonAsciiPath(unittest.TestCase):

    @onlyif_unicode_paths
    def setUp(self):
        self.BASETESTDIR = tempfile.mkdtemp()
        self.TESTDIR = join(self.BASETESTDIR, u"åäö")
        os.mkdir(self.TESTDIR)
        with open(join(self.TESTDIR, u"åäötestscript.py"), "w") as sfile:
            sfile.write("pass\n")
        self.oldpath = os.getcwd()
        os.chdir(self.TESTDIR)
        self.fname = u"åäötestscript.py"

    def tearDown(self):
        os.chdir(self.oldpath)
        shutil.rmtree(self.BASETESTDIR)

    @onlyif_unicode_paths
    def test_1(self):
        """Test safe_execfile with non-ascii path
        """
        ip.safe_execfile(self.fname, {}, raise_exceptions=True)

class ExitCodeChecks(tt.TempFileMixin):
    def test_exit_code_ok(self):
        self.system('exit 0')
        self.assertEqual(ip.user_ns['_exit_code'], 0)

    def test_exit_code_error(self):
        self.system('exit 1')
        self.assertEqual(ip.user_ns['_exit_code'], 1)
    
    @skipif(not hasattr(signal, 'SIGALRM'))
    def test_exit_code_signal(self):
        self.mktmp("import signal, time\n"
                   "signal.setitimer(signal.ITIMER_REAL, 0.1)\n"
                   "time.sleep(1)\n")
        self.system("%s %s" % (sys.executable, self.fname))
        self.assertEqual(ip.user_ns['_exit_code'], -signal.SIGALRM)
    
    @onlyif_cmds_exist("csh")
    def test_exit_code_signal_csh(self):
        SHELL = os.environ.get('SHELL', None)
        os.environ['SHELL'] = find_cmd("csh")
        try:
            self.test_exit_code_signal()
        finally:
            if SHELL is not None:
                os.environ['SHELL'] = SHELL
            else:
                del os.environ['SHELL']

class TestSystemRaw(unittest.TestCase, ExitCodeChecks):
    system = ip.system_raw

    @onlyif_unicode_paths
    def test_1(self):
        """Test system_raw with non-ascii cmd
        """
        cmd = u'''python -c "'åäö'"   '''
        ip.system_raw(cmd)

    @mock.patch('subprocess.call', side_effect=KeyboardInterrupt)
    @mock.patch('os.system', side_effect=KeyboardInterrupt)
    def test_control_c(self, *mocks):
        try:
            self.system("sleep 1 # wont happen")
        except KeyboardInterrupt:
            self.fail("system call should intercept "
                      "keyboard interrupt from subprocess.call")
        self.assertEqual(ip.user_ns['_exit_code'], -signal.SIGINT)

# TODO: Exit codes are currently ignored on Windows.
class TestSystemPipedExitCode(unittest.TestCase, ExitCodeChecks):
    system = ip.system_piped

    @skip_win32
    def test_exit_code_ok(self):
        ExitCodeChecks.test_exit_code_ok(self)

    @skip_win32
    def test_exit_code_error(self):
        ExitCodeChecks.test_exit_code_error(self)

    @skip_win32
    def test_exit_code_signal(self):
        ExitCodeChecks.test_exit_code_signal(self)

class TestModules(unittest.TestCase, tt.TempFileMixin):
    def test_extraneous_loads(self):
        """Test we're not loading modules on startup that we shouldn't.
        """
        self.mktmp("import sys\n"
                   "print('numpy' in sys.modules)\n"
                   "print('ipyparallel' in sys.modules)\n"
                   "print('ipykernel' in sys.modules)\n"
                   )
        out = "False\nFalse\nFalse\n"
        tt.ipexec_validate(self.fname, out)

class Negator(ast.NodeTransformer):
    """Negates all number literals in an AST."""
    def visit_Num(self, node):
        node.n = -node.n
        return node

class TestAstTransform(unittest.TestCase):
    def setUp(self):
        self.negator = Negator()
        ip.ast_transformers.append(self.negator)
    
    def tearDown(self):
        ip.ast_transformers.remove(self.negator)
    
    def test_run_cell(self):
        with tt.AssertPrints('-34'):
            ip.run_cell('print (12 + 22)')
        
        # A named reference to a number shouldn't be transformed.
        ip.user_ns['n'] = 55
        with tt.AssertNotPrints('-55'):
            ip.run_cell('print (n)')
    
    def test_timeit(self):
        called = set()
        def f(x):
            called.add(x)
        ip.push({'f':f})
        
        with tt.AssertPrints("std. dev. of"):
            ip.run_line_magic("timeit", "-n1 f(1)")
        self.assertEqual(called, {-1})
        called.clear()

        with tt.AssertPrints("std. dev. of"):
            ip.run_cell_magic("timeit", "-n1 f(2)", "f(3)")
        self.assertEqual(called, {-2, -3})
    
    def test_time(self):
        called = []
        def f(x):
            called.append(x)
        ip.push({'f':f})
        
        # Test with an expression
        with tt.AssertPrints("Wall time: "):
            ip.run_line_magic("time", "f(5+9)")
        self.assertEqual(called, [-14])
        called[:] = []
        
        # Test with a statement (different code path)
        with tt.AssertPrints("Wall time: "):
            ip.run_line_magic("time", "a = f(-3 + -2)")
        self.assertEqual(called, [5])
    
    def test_macro(self):
        ip.push({'a':10})
        # The AST transformation makes this do a+=-1
        ip.define_macro("amacro", "a+=1\nprint(a)")
        
        with tt.AssertPrints("9"):
            ip.run_cell("amacro")
        with tt.AssertPrints("8"):
            ip.run_cell("amacro")

class IntegerWrapper(ast.NodeTransformer):
    """Wraps all integers in a call to Integer()"""
    def visit_Num(self, node):
        if isinstance(node.n, int):
            return ast.Call(func=ast.Name(id='Integer', ctx=ast.Load()),
                            args=[node], keywords=[])
        return node

class TestAstTransform2(unittest.TestCase):
    def setUp(self):
        self.intwrapper = IntegerWrapper()
        ip.ast_transformers.append(self.intwrapper)
        
        self.calls = []
        def Integer(*args):
            self.calls.append(args)
            return args
        ip.push({"Integer": Integer})
    
    def tearDown(self):
        ip.ast_transformers.remove(self.intwrapper)
        del ip.user_ns['Integer']
    
    def test_run_cell(self):
        ip.run_cell("n = 2")
        self.assertEqual(self.calls, [(2,)])
        
        # This shouldn't throw an error
        ip.run_cell("o = 2.0")
        self.assertEqual(ip.user_ns['o'], 2.0)
    
    def test_timeit(self):
        called = set()
        def f(x):
            called.add(x)
        ip.push({'f':f})

        with tt.AssertPrints("std. dev. of"):
            ip.run_line_magic("timeit", "-n1 f(1)")
        self.assertEqual(called, {(1,)})
        called.clear()

        with tt.AssertPrints("std. dev. of"):
            ip.run_cell_magic("timeit", "-n1 f(2)", "f(3)")
        self.assertEqual(called, {(2,), (3,)})

class ErrorTransformer(ast.NodeTransformer):
    """Throws an error when it sees a number."""
    def visit_Num(self, node):
        raise ValueError("test")

class TestAstTransformError(unittest.TestCase):
    def test_unregistering(self):
        err_transformer = ErrorTransformer()
        ip.ast_transformers.append(err_transformer)
        
        with tt.AssertPrints("unregister", channel='stderr'):
            ip.run_cell("1 + 2")
        
        # This should have been removed.
        nt.assert_not_in(err_transformer, ip.ast_transformers)


class StringRejector(ast.NodeTransformer):
    """Throws an InputRejected when it sees a string literal.

    Used to verify that NodeTransformers can signal that a piece of code should
    not be executed by throwing an InputRejected.
    """

    def visit_Str(self, node):
        raise InputRejected("test")


class TestAstTransformInputRejection(unittest.TestCase):

    def setUp(self):
        self.transformer = StringRejector()
        ip.ast_transformers.append(self.transformer)

    def tearDown(self):
        ip.ast_transformers.remove(self.transformer)

    def test_input_rejection(self):
        """Check that NodeTransformers can reject input."""

        expect_exception_tb = tt.AssertPrints("InputRejected: test")
        expect_no_cell_output = tt.AssertNotPrints("'unsafe'", suppress=False)

        # Run the same check twice to verify that the transformer is not
        # disabled after raising.
        with expect_exception_tb, expect_no_cell_output:
            ip.run_cell("'unsafe'")

        with expect_exception_tb, expect_no_cell_output:
            res = ip.run_cell("'unsafe'")

        self.assertIsInstance(res.error_before_exec, InputRejected)

def test__IPYTHON__():
    # This shouldn't raise a NameError, that's all
    __IPYTHON__


class DummyRepr(object):
    def __repr__(self):
        return "DummyRepr"
    
    def _repr_html_(self):
        return "<b>dummy</b>"
    
    def _repr_javascript_(self):
        return "console.log('hi');", {'key': 'value'}
    

def test_user_variables():
    # enable all formatters
    ip.display_formatter.active_types = ip.display_formatter.format_types
    
    ip.user_ns['dummy'] = d = DummyRepr()
    keys = {'dummy', 'doesnotexist'}
    r = ip.user_expressions({ key:key for key in keys})

    nt.assert_equal(keys, set(r.keys()))
    dummy = r['dummy']
    nt.assert_equal({'status', 'data', 'metadata'}, set(dummy.keys()))
    nt.assert_equal(dummy['status'], 'ok')
    data = dummy['data']
    metadata = dummy['metadata']
    nt.assert_equal(data.get('text/html'), d._repr_html_())
    js, jsmd = d._repr_javascript_()
    nt.assert_equal(data.get('application/javascript'), js)
    nt.assert_equal(metadata.get('application/javascript'), jsmd)
    
    dne = r['doesnotexist']
    nt.assert_equal(dne['status'], 'error')
    nt.assert_equal(dne['ename'], 'NameError')
    
    # back to text only
    ip.display_formatter.active_types = ['text/plain']
    
def test_user_expression():
    # enable all formatters
    ip.display_formatter.active_types = ip.display_formatter.format_types
    query = {
        'a' : '1 + 2',
        'b' : '1/0',
    }
    r = ip.user_expressions(query)
    import pprint
    pprint.pprint(r)
    nt.assert_equal(set(r.keys()), set(query.keys()))
    a = r['a']
    nt.assert_equal({'status', 'data', 'metadata'}, set(a.keys()))
    nt.assert_equal(a['status'], 'ok')
    data = a['data']
    metadata = a['metadata']
    nt.assert_equal(data.get('text/plain'), '3')
    
    b = r['b']
    nt.assert_equal(b['status'], 'error')
    nt.assert_equal(b['ename'], 'ZeroDivisionError')
    
    # back to text only
    ip.display_formatter.active_types = ['text/plain']
    




class TestSyntaxErrorTransformer(unittest.TestCase):
    """Check that SyntaxError raised by an input transformer is handled by run_cell()"""

    class SyntaxErrorTransformer(InputTransformer):

        def push(self, line):
            pos = line.find('syntaxerror')
            if pos >= 0:
                e = SyntaxError('input contains "syntaxerror"')
                e.text = line
                e.offset = pos + 1
                raise e
            return line

        def reset(self):
            pass

    def setUp(self):
        self.transformer = TestSyntaxErrorTransformer.SyntaxErrorTransformer()
        ip.input_splitter.python_line_transforms.append(self.transformer)
        ip.input_transformer_manager.python_line_transforms.append(self.transformer)

    def tearDown(self):
        ip.input_splitter.python_line_transforms.remove(self.transformer)
        ip.input_transformer_manager.python_line_transforms.remove(self.transformer)

    def test_syntaxerror_input_transformer(self):
        with tt.AssertPrints('1234'):
            ip.run_cell('1234')
        with tt.AssertPrints('SyntaxError: invalid syntax'):
            ip.run_cell('1 2 3')   # plain python syntax error
        with tt.AssertPrints('SyntaxError: input contains "syntaxerror"'):
            ip.run_cell('2345  # syntaxerror')  # input transformer syntax error
        with tt.AssertPrints('3456'):
            ip.run_cell('3456')



def test_warning_suppression():
    ip.run_cell("import warnings")
    try:
        with tt.AssertPrints("UserWarning: asdf", channel="stderr"):
            ip.run_cell("warnings.warn('asdf')")
        # Here's the real test -- if we run that again, we should get the
        # warning again. Traditionally, each warning was only issued once per
        # IPython session (approximately), even if the user typed in new and
        # different code that should have also triggered the warning, leading
        # to much confusion.
        with tt.AssertPrints("UserWarning: asdf", channel="stderr"):
            ip.run_cell("warnings.warn('asdf')")
    finally:
        ip.run_cell("del warnings")


def test_deprecation_warning():
    ip.run_cell("""
import warnings
def wrn():
    warnings.warn(
        "I AM  A WARNING",
        DeprecationWarning
    )
        """)
    try:
        with tt.AssertPrints("I AM  A WARNING", channel="stderr"):
            ip.run_cell("wrn()")
    finally:
        ip.run_cell("del warnings")
        ip.run_cell("del wrn")


class TestImportNoDeprecate(tt.TempFileMixin):

    def setup(self):
        """Make a valid python temp file."""
        self.mktmp("""
import warnings
def wrn():
    warnings.warn(
        "I AM  A WARNING",
        DeprecationWarning
    )
""")

    def test_no_dep(self):
        """
        No deprecation warning should be raised from imported functions
        """
        ip.run_cell("from {} import wrn".format(self.fname))

        with tt.AssertNotPrints("I AM  A WARNING"):
            ip.run_cell("wrn()")
        ip.run_cell("del wrn")


def test_custom_exc_count():
    hook = mock.Mock(return_value=None)
    ip.set_custom_exc((SyntaxError,), hook)
    before = ip.execution_count
    ip.run_cell("def foo()", store_history=True)
    # restore default excepthook
    ip.set_custom_exc((), None)
    nt.assert_equal(hook.call_count, 1)
    nt.assert_equal(ip.execution_count, before + 1)
