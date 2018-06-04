# -*- coding: utf-8 -*-
"""Tests for the TerminalInteractiveShell and related pieces."""
# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import sys
import unittest

from IPython.core.inputtransformer import InputTransformer
from IPython.testing import tools as tt
from IPython.utils.capture import capture_output

from IPython.terminal.ptutils import _elide, _adjust_completion_text_based_on_context
import nose.tools as nt

class TestElide(unittest.TestCase):

    def test_elide(self):
        _elide('concatenate((a1, a2, ...), axis') # do not raise
        _elide('concatenate((a1, a2, ..), . axis') # do not raise
        nt.assert_equal(_elide('aaaa.bbbb.ccccc.dddddd.eeeee.fffff.gggggg.hhhhhh'), 'aaaa.b…g.hhhhhh')


class TestContextAwareCompletion(unittest.TestCase):

    def test_adjust_completion_text_based_on_context(self):
        # Adjusted case
        nt.assert_equal(_adjust_completion_text_based_on_context('arg1=', 'func1(a=)', 7), 'arg1')

        # Untouched cases
        nt.assert_equal(_adjust_completion_text_based_on_context('arg1=', 'func1(a)', 7), 'arg1=')
        nt.assert_equal(_adjust_completion_text_based_on_context('arg1=', 'func1(a', 7), 'arg1=')
        nt.assert_equal(_adjust_completion_text_based_on_context('%magic', 'func1(a=)', 7), '%magic')
        nt.assert_equal(_adjust_completion_text_based_on_context('func2', 'func1(a=)', 7), 'func2')

# Decorator for interaction loop tests -----------------------------------------

class mock_input_helper(object):
    """Machinery for tests of the main interact loop.

    Used by the mock_input decorator.
    """
    def __init__(self, testgen):
        self.testgen = testgen
        self.exception = None
        self.ip = get_ipython()

    def __enter__(self):
        self.orig_prompt_for_code = self.ip.prompt_for_code
        self.ip.prompt_for_code = self.fake_input
        return self

    def __exit__(self, etype, value, tb):
        self.ip.prompt_for_code = self.orig_prompt_for_code

    def fake_input(self):
        try:
            return next(self.testgen)
        except StopIteration:
            self.ip.keep_running = False
            return u''
        except:
            self.exception = sys.exc_info()
            self.ip.keep_running = False
            return u''

def mock_input(testfunc):
    """Decorator for tests of the main interact loop.

    Write the test as a generator, yield-ing the input strings, which IPython
    will see as if they were typed in at the prompt.
    """
    def test_method(self):
        testgen = testfunc(self)
        with mock_input_helper(testgen) as mih:
            mih.ip.interact()

        if mih.exception is not None:
            # Re-raise captured exception
            etype, value, tb = mih.exception
            import traceback
            traceback.print_tb(tb, file=sys.stdout)
            del tb  # Avoid reference loop
            raise value

    return test_method

# Test classes -----------------------------------------------------------------

class InteractiveShellTestCase(unittest.TestCase):
    def rl_hist_entries(self, rl, n):
        """Get last n readline history entries as a list"""
        return [rl.get_history_item(rl.get_current_history_length() - x)
                for x in range(n - 1, -1, -1)]
    
    @mock_input
    def test_inputtransformer_syntaxerror(self):
        ip = get_ipython()
        transformer = SyntaxErrorTransformer()
        ip.input_splitter.python_line_transforms.append(transformer)
        ip.input_transformer_manager.python_line_transforms.append(transformer)

        try:
            #raise Exception
            with tt.AssertPrints('4', suppress=False):
                yield u'print(2*2)'

            with tt.AssertPrints('SyntaxError: input contains', suppress=False):
                yield u'print(2345) # syntaxerror'

            with tt.AssertPrints('16', suppress=False):
                yield u'print(4*4)'

        finally:
            ip.input_splitter.python_line_transforms.remove(transformer)
            ip.input_transformer_manager.python_line_transforms.remove(transformer)

    def test_plain_text_only(self):
        ip = get_ipython()
        formatter = ip.display_formatter
        assert formatter.active_types == ['text/plain']
        assert not formatter.ipython_display_formatter.enabled

        class Test(object):
            def __repr__(self):
                return "<Test %i>" % id(self)

            def _repr_html_(self):
                return '<html>'

        # verify that HTML repr isn't computed
        obj = Test()
        data, _ = formatter.format(obj)
        self.assertEqual(data, {'text/plain': repr(obj)})

        class Test2(Test):
            def _ipython_display_(self):
                from IPython.display import display
                display('<custom>')

        # verify that _ipython_display_ shortcut isn't called
        obj = Test2()
        with capture_output() as captured:
            data, _ = formatter.format(obj)

        self.assertEqual(data, {'text/plain': repr(obj)})
        assert captured.stdout == ''



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

class TerminalMagicsTestCase(unittest.TestCase):
    def test_paste_magics_blankline(self):
        """Test that code with a blank line doesn't get split (gh-3246)."""
        ip = get_ipython()
        s = ('def pasted_func(a):\n'
             '    b = a+1\n'
             '\n'
             '    return b')
        
        tm = ip.magics_manager.registry['TerminalMagics']
        tm.store_or_execute(s, name=None)
        
        self.assertEqual(ip.user_ns['pasted_func'](54), 55)
