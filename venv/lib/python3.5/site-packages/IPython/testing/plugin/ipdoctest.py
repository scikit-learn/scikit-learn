"""Nose Plugin that supports IPython doctests.

Limitations:

- When generating examples for use as doctests, make sure that you have
  pretty-printing OFF.  This can be done either by setting the
  ``PlainTextFormatter.pprint`` option in your configuration file to  False, or
  by interactively disabling it with  %Pprint.  This is required so that IPython
  output matches that of normal Python, which is used by doctest for internal
  execution.

- Do not rely on specific prompt numbers for results (such as using
  '_34==True', for example).  For IPython tests run via an external process the
  prompt numbers may be different, and IPython tests run as normal python code
  won't even have these special _NN variables set at all.
"""

#-----------------------------------------------------------------------------
# Module imports

# From the standard library
import builtins as builtin_mod
import doctest
import inspect
import logging
import os
import re
import sys
from importlib import import_module
from io import StringIO

from testpath import modified_env

from inspect import getmodule

# We are overriding the default doctest runner, so we need to import a few
# things from doctest directly
from doctest import (REPORTING_FLAGS, REPORT_ONLY_FIRST_FAILURE,
                     _unittest_reportflags, DocTestRunner,
                     _extract_future_flags, pdb, _OutputRedirectingPdb,
                     _exception_traceback,
                     linecache)

# Third-party modules

from nose.plugins import doctests, Plugin
from nose.util import anyp, tolist

#-----------------------------------------------------------------------------
# Module globals and other constants
#-----------------------------------------------------------------------------

log = logging.getLogger(__name__)


#-----------------------------------------------------------------------------
# Classes and functions
#-----------------------------------------------------------------------------

def is_extension_module(filename):
    """Return whether the given filename is an extension module.

    This simply checks that the extension is either .so or .pyd.
    """
    return os.path.splitext(filename)[1].lower() in ('.so','.pyd')


class DocTestSkip(object):
    """Object wrapper for doctests to be skipped."""

    ds_skip = """Doctest to skip.
    >>> 1 #doctest: +SKIP
    """

    def __init__(self,obj):
        self.obj = obj

    def __getattribute__(self,key):
        if key == '__doc__':
            return DocTestSkip.ds_skip
        else:
            return getattr(object.__getattribute__(self,'obj'),key)

# Modified version of the one in the stdlib, that fixes a python bug (doctests
# not found in extension modules, http://bugs.python.org/issue3158)
class DocTestFinder(doctest.DocTestFinder):

    def _from_module(self, module, object):
        """
        Return true if the given object is defined in the given
        module.
        """
        if module is None:
            return True
        elif inspect.isfunction(object):
            return module.__dict__ is object.__globals__
        elif inspect.isbuiltin(object):
            return module.__name__ == object.__module__
        elif inspect.isclass(object):
            return module.__name__ == object.__module__
        elif inspect.ismethod(object):
            # This one may be a bug in cython that fails to correctly set the
            # __module__ attribute of methods, but since the same error is easy
            # to make by extension code writers, having this safety in place
            # isn't such a bad idea
            return module.__name__ == object.__self__.__class__.__module__
        elif inspect.getmodule(object) is not None:
            return module is inspect.getmodule(object)
        elif hasattr(object, '__module__'):
            return module.__name__ == object.__module__
        elif isinstance(object, property):
            return True # [XX] no way not be sure.
        elif inspect.ismethoddescriptor(object):
            # Unbound PyQt signals reach this point in Python 3.4b3, and we want
            # to avoid throwing an error. See also http://bugs.python.org/issue3158
            return False
        else:
            raise ValueError("object must be a class or function, got %r" % object)

    def _find(self, tests, obj, name, module, source_lines, globs, seen):
        """
        Find tests for the given object and any contained objects, and
        add them to `tests`.
        """
        print('_find for:', obj, name, module)  # dbg
        if hasattr(obj,"skip_doctest"):
            #print 'SKIPPING DOCTEST FOR:',obj  # dbg
            obj = DocTestSkip(obj)
        
        doctest.DocTestFinder._find(self,tests, obj, name, module,
                                    source_lines, globs, seen)

        # Below we re-run pieces of the above method with manual modifications,
        # because the original code is buggy and fails to correctly identify
        # doctests in extension modules.

        # Local shorthands
        from inspect import isroutine, isclass

        # Look for tests in a module's contained objects.
        if inspect.ismodule(obj) and self._recurse:
            for valname, val in obj.__dict__.items():
                valname1 = '%s.%s' % (name, valname)
                if ( (isroutine(val) or isclass(val))
                     and self._from_module(module, val) ):

                    self._find(tests, val, valname1, module, source_lines,
                               globs, seen)

        # Look for tests in a class's contained objects.
        if inspect.isclass(obj) and self._recurse:
            #print 'RECURSE into class:',obj  # dbg
            for valname, val in obj.__dict__.items():
                # Special handling for staticmethod/classmethod.
                if isinstance(val, staticmethod):
                    val = getattr(obj, valname)
                if isinstance(val, classmethod):
                    val = getattr(obj, valname).__func__

                # Recurse to methods, properties, and nested classes.
                if ((inspect.isfunction(val) or inspect.isclass(val) or
                     inspect.ismethod(val) or
                      isinstance(val, property)) and
                      self._from_module(module, val)):
                    valname = '%s.%s' % (name, valname)
                    self._find(tests, val, valname, module, source_lines,
                               globs, seen)


class IPDoctestOutputChecker(doctest.OutputChecker):
    """Second-chance checker with support for random tests.

    If the default comparison doesn't pass, this checker looks in the expected
    output string for flags that tell us to ignore the output.
    """

    random_re = re.compile(r'#\s*random\s+')

    def check_output(self, want, got, optionflags):
        """Check output, accepting special markers embedded in the output.

        If the output didn't pass the default validation but the special string
        '#random' is included, we accept it."""

        # Let the original tester verify first, in case people have valid tests
        # that happen to have a comment saying '#random' embedded in.
        ret = doctest.OutputChecker.check_output(self, want, got,
                                                 optionflags)
        if not ret and self.random_re.search(want):
            #print >> sys.stderr, 'RANDOM OK:',want  # dbg
            return True

        return ret


class DocTestCase(doctests.DocTestCase):
    """Proxy for DocTestCase: provides an address() method that
    returns the correct address for the doctest case. Otherwise
    acts as a proxy to the test case. To provide hints for address(),
    an obj may also be passed -- this will be used as the test object
    for purposes of determining the test address, if it is provided.
    """

    # Note: this method was taken from numpy's nosetester module.

    # Subclass nose.plugins.doctests.DocTestCase to work around a bug in
    # its constructor that blocks non-default arguments from being passed
    # down into doctest.DocTestCase

    def __init__(self, test, optionflags=0, setUp=None, tearDown=None,
                 checker=None, obj=None, result_var='_'):
        self._result_var = result_var
        doctests.DocTestCase.__init__(self, test,
                                      optionflags=optionflags,
                                      setUp=setUp, tearDown=tearDown,
                                      checker=checker)
        # Now we must actually copy the original constructor from the stdlib
        # doctest class, because we can't call it directly and a bug in nose
        # means it never gets passed the right arguments.

        self._dt_optionflags = optionflags
        self._dt_checker = checker
        self._dt_test = test
        self._dt_test_globs_ori = test.globs
        self._dt_setUp = setUp
        self._dt_tearDown = tearDown

        # XXX - store this runner once in the object!
        runner = IPDocTestRunner(optionflags=optionflags,
                                 checker=checker, verbose=False)
        self._dt_runner = runner


        # Each doctest should remember the directory it was loaded from, so
        # things like %run work without too many contortions
        self._ori_dir = os.path.dirname(test.filename)

    # Modified runTest from the default stdlib
    def runTest(self):
        test = self._dt_test
        runner = self._dt_runner

        old = sys.stdout
        new = StringIO()
        optionflags = self._dt_optionflags

        if not (optionflags & REPORTING_FLAGS):
            # The option flags don't include any reporting flags,
            # so add the default reporting flags
            optionflags |= _unittest_reportflags

        try:
            # Save our current directory and switch out to the one where the
            # test was originally created, in case another doctest did a
            # directory change.  We'll restore this in the finally clause.
            curdir = os.getcwd()
            #print 'runTest in dir:', self._ori_dir  # dbg
            os.chdir(self._ori_dir)

            runner.DIVIDER = "-"*70
            failures, tries = runner.run(test,out=new.write,
                                         clear_globs=False)
        finally:
            sys.stdout = old
            os.chdir(curdir)

        if failures:
            raise self.failureException(self.format_failure(new.getvalue()))

    def setUp(self):
        """Modified test setup that syncs with ipython namespace"""
        #print "setUp test", self._dt_test.examples # dbg
        if isinstance(self._dt_test.examples[0], IPExample):
            # for IPython examples *only*, we swap the globals with the ipython
            # namespace, after updating it with the globals (which doctest
            # fills with the necessary info from the module being tested).
            self.user_ns_orig = {}
            self.user_ns_orig.update(_ip.user_ns)
            _ip.user_ns.update(self._dt_test.globs)
            # We must remove the _ key in the namespace, so that Python's
            # doctest code sets it naturally
            _ip.user_ns.pop('_', None)
            _ip.user_ns['__builtins__'] = builtin_mod
            self._dt_test.globs = _ip.user_ns

        super(DocTestCase, self).setUp()

    def tearDown(self):

        # Undo the test.globs reassignment we made, so that the parent class
        # teardown doesn't destroy the ipython namespace
        if isinstance(self._dt_test.examples[0], IPExample):
            self._dt_test.globs = self._dt_test_globs_ori
            _ip.user_ns.clear()
            _ip.user_ns.update(self.user_ns_orig)

        # XXX - fperez: I am not sure if this is truly a bug in nose 0.11, but
        # it does look like one to me: its tearDown method tries to run
        #
        # delattr(builtin_mod, self._result_var)
        #
        # without checking that the attribute really is there; it implicitly
        # assumes it should have been set via displayhook.  But if the
        # displayhook was never called, this doesn't necessarily happen.  I
        # haven't been able to find a little self-contained example outside of
        # ipython that would show the problem so I can report it to the nose
        # team, but it does happen a lot in our code.
        #
        # So here, we just protect as narrowly as possible by trapping an
        # attribute error whose message would be the name of self._result_var,
        # and letting any other error propagate.
        try:
            super(DocTestCase, self).tearDown()
        except AttributeError as exc:
            if exc.args[0] != self._result_var:
                raise


# A simple subclassing of the original with a different class name, so we can
# distinguish and treat differently IPython examples from pure python ones.
class IPExample(doctest.Example): pass


class IPExternalExample(doctest.Example):
    """Doctest examples to be run in an external process."""

    def __init__(self, source, want, exc_msg=None, lineno=0, indent=0,
                 options=None):
        # Parent constructor
        doctest.Example.__init__(self,source,want,exc_msg,lineno,indent,options)

        # An EXTRA newline is needed to prevent pexpect hangs
        self.source += '\n'


class IPDocTestParser(doctest.DocTestParser):
    """
    A class used to parse strings containing doctest examples.

    Note: This is a version modified to properly recognize IPython input and
    convert any IPython examples into valid Python ones.
    """
    # This regular expression is used to find doctest examples in a
    # string.  It defines three groups: `source` is the source code
    # (including leading indentation and prompts); `indent` is the
    # indentation of the first (PS1) line of the source code; and
    # `want` is the expected output (including leading indentation).

    # Classic Python prompts or default IPython ones
    _PS1_PY = r'>>>'
    _PS2_PY = r'\.\.\.'

    _PS1_IP = r'In\ \[\d+\]:'
    _PS2_IP = r'\ \ \ \.\.\.+:'

    _RE_TPL = r'''
        # Source consists of a PS1 line followed by zero or more PS2 lines.
        (?P<source>
            (?:^(?P<indent> [ ]*) (?P<ps1> %s) .*)    # PS1 line
            (?:\n           [ ]*  (?P<ps2> %s) .*)*)  # PS2 lines
        \n? # a newline
        # Want consists of any non-blank lines that do not start with PS1.
        (?P<want> (?:(?![ ]*$)    # Not a blank line
                     (?![ ]*%s)   # Not a line starting with PS1
                     (?![ ]*%s)   # Not a line starting with PS2
                     .*$\n?       # But any other line
                  )*)
                  '''

    _EXAMPLE_RE_PY = re.compile( _RE_TPL % (_PS1_PY,_PS2_PY,_PS1_PY,_PS2_PY),
                                 re.MULTILINE | re.VERBOSE)

    _EXAMPLE_RE_IP = re.compile( _RE_TPL % (_PS1_IP,_PS2_IP,_PS1_IP,_PS2_IP),
                                 re.MULTILINE | re.VERBOSE)

    # Mark a test as being fully random.  In this case, we simply append the
    # random marker ('#random') to each individual example's output.  This way
    # we don't need to modify any other code.
    _RANDOM_TEST = re.compile(r'#\s*all-random\s+')

    # Mark tests to be executed in an external process - currently unsupported.
    _EXTERNAL_IP = re.compile(r'#\s*ipdoctest:\s*EXTERNAL')

    def ip2py(self,source):
        """Convert input IPython source into valid Python."""
        block = _ip.input_transformer_manager.transform_cell(source)
        if len(block.splitlines()) == 1:
            return _ip.prefilter(block)
        else:
            return block

    def parse(self, string, name='<string>'):
        """
        Divide the given string into examples and intervening text,
        and return them as a list of alternating Examples and strings.
        Line numbers for the Examples are 0-based.  The optional
        argument `name` is a name identifying this string, and is only
        used for error messages.
        """

        #print 'Parse string:\n',string # dbg

        string = string.expandtabs()
        # If all lines begin with the same indentation, then strip it.
        min_indent = self._min_indent(string)
        if min_indent > 0:
            string = '\n'.join([l[min_indent:] for l in string.split('\n')])

        output = []
        charno, lineno = 0, 0

        # We make 'all random' tests by adding the '# random' mark to every
        # block of output in the test.
        if self._RANDOM_TEST.search(string):
            random_marker = '\n# random'
        else:
            random_marker = ''

        # Whether to convert the input from ipython to python syntax
        ip2py = False
        # Find all doctest examples in the string.  First, try them as Python
        # examples, then as IPython ones
        terms = list(self._EXAMPLE_RE_PY.finditer(string))
        if terms:
            # Normal Python example
            #print '-'*70  # dbg
            #print 'PyExample, Source:\n',string  # dbg
            #print '-'*70  # dbg
            Example = doctest.Example
        else:
            # It's an ipython example.  Note that IPExamples are run
            # in-process, so their syntax must be turned into valid python.
            # IPExternalExamples are run out-of-process (via pexpect) so they
            # don't need any filtering (a real ipython will be executing them).
            terms = list(self._EXAMPLE_RE_IP.finditer(string))
            if self._EXTERNAL_IP.search(string):
                #print '-'*70  # dbg
                #print 'IPExternalExample, Source:\n',string  # dbg
                #print '-'*70  # dbg
                Example = IPExternalExample
            else:
                #print '-'*70  # dbg
                #print 'IPExample, Source:\n',string  # dbg
                #print '-'*70  # dbg
                Example = IPExample
                ip2py = True

        for m in terms:
            # Add the pre-example text to `output`.
            output.append(string[charno:m.start()])
            # Update lineno (lines before this example)
            lineno += string.count('\n', charno, m.start())
            # Extract info from the regexp match.
            (source, options, want, exc_msg) = \
                     self._parse_example(m, name, lineno,ip2py)

            # Append the random-output marker (it defaults to empty in most
            # cases, it's only non-empty for 'all-random' tests):
            want += random_marker

            if Example is IPExternalExample:
                options[doctest.NORMALIZE_WHITESPACE] = True
                want += '\n'

            # Create an Example, and add it to the list.
            if not self._IS_BLANK_OR_COMMENT(source):
                output.append(Example(source, want, exc_msg,
                                      lineno=lineno,
                                      indent=min_indent+len(m.group('indent')),
                                      options=options))
            # Update lineno (lines inside this example)
            lineno += string.count('\n', m.start(), m.end())
            # Update charno.
            charno = m.end()
        # Add any remaining post-example text to `output`.
        output.append(string[charno:])
        return output

    def _parse_example(self, m, name, lineno,ip2py=False):
        """
        Given a regular expression match from `_EXAMPLE_RE` (`m`),
        return a pair `(source, want)`, where `source` is the matched
        example's source code (with prompts and indentation stripped);
        and `want` is the example's expected output (with indentation
        stripped).

        `name` is the string's name, and `lineno` is the line number
        where the example starts; both are used for error messages.

        Optional:
        `ip2py`: if true, filter the input via IPython to convert the syntax
        into valid python.
        """

        # Get the example's indentation level.
        indent = len(m.group('indent'))

        # Divide source into lines; check that they're properly
        # indented; and then strip their indentation & prompts.
        source_lines = m.group('source').split('\n')

        # We're using variable-length input prompts
        ps1 = m.group('ps1')
        ps2 = m.group('ps2')
        ps1_len = len(ps1)

        self._check_prompt_blank(source_lines, indent, name, lineno,ps1_len)
        if ps2:
            self._check_prefix(source_lines[1:], ' '*indent + ps2, name, lineno)

        source = '\n'.join([sl[indent+ps1_len+1:] for sl in source_lines])

        if ip2py:
            # Convert source input from IPython into valid Python syntax
            source = self.ip2py(source)

        # Divide want into lines; check that it's properly indented; and
        # then strip the indentation.  Spaces before the last newline should
        # be preserved, so plain rstrip() isn't good enough.
        want = m.group('want')
        want_lines = want.split('\n')
        if len(want_lines) > 1 and re.match(r' *$', want_lines[-1]):
            del want_lines[-1]  # forget final newline & spaces after it
        self._check_prefix(want_lines, ' '*indent, name,
                           lineno + len(source_lines))

        # Remove ipython output prompt that might be present in the first line
        want_lines[0] = re.sub(r'Out\[\d+\]: \s*?\n?','',want_lines[0])

        want = '\n'.join([wl[indent:] for wl in want_lines])

        # If `want` contains a traceback message, then extract it.
        m = self._EXCEPTION_RE.match(want)
        if m:
            exc_msg = m.group('msg')
        else:
            exc_msg = None

        # Extract options from the source.
        options = self._find_options(source, name, lineno)

        return source, options, want, exc_msg

    def _check_prompt_blank(self, lines, indent, name, lineno, ps1_len):
        """
        Given the lines of a source string (including prompts and
        leading indentation), check to make sure that every prompt is
        followed by a space character.  If any line is not followed by
        a space character, then raise ValueError.

        Note: IPython-modified version which takes the input prompt length as a
        parameter, so that prompts of variable length can be dealt with.
        """
        space_idx = indent+ps1_len
        min_len = space_idx+1
        for i, line in enumerate(lines):
            if len(line) >=  min_len and line[space_idx] != ' ':
                raise ValueError('line %r of the docstring for %s '
                                 'lacks blank after %s: %r' %
                                 (lineno+i+1, name,
                                  line[indent:space_idx], line))


SKIP = doctest.register_optionflag('SKIP')


class IPDocTestRunner(doctest.DocTestRunner,object):
    """Test runner that synchronizes the IPython namespace with test globals.
    """

    def run(self, test, compileflags=None, out=None, clear_globs=True):

        # Hack: ipython needs access to the execution context of the example,
        # so that it can propagate user variables loaded by %run into
        # test.globs.  We put them here into our modified %run as a function
        # attribute.  Our new %run will then only make the namespace update
        # when called (rather than unconconditionally updating test.globs here
        # for all examples, most of which won't be calling %run anyway).
        #_ip._ipdoctest_test_globs = test.globs
        #_ip._ipdoctest_test_filename = test.filename

        test.globs.update(_ip.user_ns)

        # Override terminal size to standardise traceback format
        with modified_env({'COLUMNS': '80', 'LINES': '24'}):
            return super(IPDocTestRunner,self).run(test,
                                                   compileflags,out,clear_globs)


class DocFileCase(doctest.DocFileCase):
    """Overrides to provide filename
    """
    def address(self):
        return (self._dt_test.filename, None, None)


class ExtensionDoctest(doctests.Doctest):
    """Nose Plugin that supports doctests in extension modules.
    """
    name = 'extdoctest'   # call nosetests with --with-extdoctest
    enabled = True

    def options(self, parser, env=os.environ):
        Plugin.options(self, parser, env)
        parser.add_option('--doctest-tests', action='store_true',
                          dest='doctest_tests',
                          default=env.get('NOSE_DOCTEST_TESTS',True),
                          help="Also look for doctests in test modules. "
                          "Note that classes, methods and functions should "
                          "have either doctests or non-doctest tests, "
                          "not both. [NOSE_DOCTEST_TESTS]")
        parser.add_option('--doctest-extension', action="append",
                          dest="doctestExtension",
                          help="Also look for doctests in files with "
                          "this extension [NOSE_DOCTEST_EXTENSION]")
        # Set the default as a list, if given in env; otherwise
        # an additional value set on the command line will cause
        # an error.
        env_setting = env.get('NOSE_DOCTEST_EXTENSION')
        if env_setting is not None:
            parser.set_defaults(doctestExtension=tolist(env_setting))


    def configure(self, options, config):
        Plugin.configure(self, options, config)
        # Pull standard doctest plugin out of config; we will do doctesting
        config.plugins.plugins = [p for p in config.plugins.plugins
                                  if p.name != 'doctest']
        self.doctest_tests = options.doctest_tests
        self.extension = tolist(options.doctestExtension)

        self.parser = doctest.DocTestParser()
        self.finder = DocTestFinder()
        self.checker = IPDoctestOutputChecker()
        self.globs = None
        self.extraglobs = None


    def loadTestsFromExtensionModule(self,filename):
        bpath,mod = os.path.split(filename)
        modname = os.path.splitext(mod)[0]
        try:
            sys.path.append(bpath)
            module = import_module(modname)
            tests = list(self.loadTestsFromModule(module))
        finally:
            sys.path.pop()
        return tests

    # NOTE: the method below is almost a copy of the original one in nose, with
    # a  few modifications to control output checking.

    def loadTestsFromModule(self, module):
        #print '*** ipdoctest - lTM',module  # dbg

        if not self.matches(module.__name__):
            log.debug("Doctest doesn't want module %s", module)
            return

        tests = self.finder.find(module,globs=self.globs,
                                 extraglobs=self.extraglobs)
        if not tests:
            return

        # always use whitespace and ellipsis options
        optionflags = doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS

        tests.sort()
        module_file = module.__file__
        if module_file[-4:] in ('.pyc', '.pyo'):
            module_file = module_file[:-1]
        for test in tests:
            if not test.examples:
                continue
            if not test.filename:
                test.filename = module_file

            yield DocTestCase(test,
                              optionflags=optionflags,
                              checker=self.checker)


    def loadTestsFromFile(self, filename):
        #print "ipdoctest - from file", filename # dbg
        if is_extension_module(filename):
            for t in self.loadTestsFromExtensionModule(filename):
                yield t
        else:
            if self.extension and anyp(filename.endswith, self.extension):
                name = os.path.basename(filename)
                dh = open(filename)
                try:
                    doc = dh.read()
                finally:
                    dh.close()
                test = self.parser.get_doctest(
                    doc, globs={'__file__': filename}, name=name,
                    filename=filename, lineno=0)
                if test.examples:
                    #print 'FileCase:',test.examples  # dbg
                    yield DocFileCase(test)
                else:
                    yield False # no tests to load


class IPythonDoctest(ExtensionDoctest):
    """Nose Plugin that supports doctests in extension modules.
    """
    name = 'ipdoctest'   # call nosetests with --with-ipdoctest
    enabled = True

    def makeTest(self, obj, parent):
        """Look for doctests in the given object, which will be a
        function, method or class.
        """
        #print 'Plugin analyzing:', obj, parent  # dbg
        # always use whitespace and ellipsis options
        optionflags = doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS

        doctests = self.finder.find(obj, module=getmodule(parent))
        if doctests:
            for test in doctests:
                if len(test.examples) == 0:
                    continue

                yield DocTestCase(test, obj=obj,
                                  optionflags=optionflags,
                                  checker=self.checker)

    def options(self, parser, env=os.environ):
        #print "Options for nose plugin:", self.name # dbg
        Plugin.options(self, parser, env)
        parser.add_option('--ipdoctest-tests', action='store_true',
                          dest='ipdoctest_tests',
                          default=env.get('NOSE_IPDOCTEST_TESTS',True),
                          help="Also look for doctests in test modules. "
                          "Note that classes, methods and functions should "
                          "have either doctests or non-doctest tests, "
                          "not both. [NOSE_IPDOCTEST_TESTS]")
        parser.add_option('--ipdoctest-extension', action="append",
                          dest="ipdoctest_extension",
                          help="Also look for doctests in files with "
                          "this extension [NOSE_IPDOCTEST_EXTENSION]")
        # Set the default as a list, if given in env; otherwise
        # an additional value set on the command line will cause
        # an error.
        env_setting = env.get('NOSE_IPDOCTEST_EXTENSION')
        if env_setting is not None:
            parser.set_defaults(ipdoctest_extension=tolist(env_setting))

    def configure(self, options, config):
        #print "Configuring nose plugin:", self.name # dbg
        Plugin.configure(self, options, config)
        # Pull standard doctest plugin out of config; we will do doctesting
        config.plugins.plugins = [p for p in config.plugins.plugins
                                  if p.name != 'doctest']
        self.doctest_tests = options.ipdoctest_tests
        self.extension = tolist(options.ipdoctest_extension)

        self.parser = IPDocTestParser()
        self.finder = DocTestFinder(parser=self.parser)
        self.checker = IPDoctestOutputChecker()
        self.globs = None
        self.extraglobs = None
