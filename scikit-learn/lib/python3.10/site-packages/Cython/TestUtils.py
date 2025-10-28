import os
import re
import unittest
import shlex
import sys
import tempfile
import textwrap
from functools import partial

from .Compiler import Errors
from .CodeWriter import CodeWriter
from .Compiler.TreeFragment import TreeFragment, strip_common_indent, StringParseContext
from .Compiler.Visitor import TreeVisitor, VisitorTransform
from .Compiler import TreePath
from .Compiler.ParseTreeTransforms import PostParse


class NodeTypeWriter(TreeVisitor):
    def __init__(self):
        super().__init__()
        self._indents = 0
        self.result = []

    def visit_Node(self, node):
        if not self.access_path:
            name = "(root)"
        else:
            tip = self.access_path[-1]
            if tip[2] is not None:
                name = "%s[%d]" % tip[1:3]
            else:
                name = tip[1]

        self.result.append("  " * self._indents +
                           "%s: %s" % (name, node.__class__.__name__))
        self._indents += 1
        self.visitchildren(node)
        self._indents -= 1


def treetypes(root):
    """Returns a string representing the tree by class names.
    There's a leading and trailing whitespace so that it can be
    compared by simple string comparison while still making test
    cases look ok."""
    w = NodeTypeWriter()
    w.visit(root)
    return "\n".join([""] + w.result + [""])


class CythonTest(unittest.TestCase):

    def setUp(self):
        Errors.init_thread()

    def tearDown(self):
        Errors.init_thread()

    def assertLines(self, expected, result):
        "Checks that the given strings or lists of strings are equal line by line"
        if not isinstance(expected, list):
            expected = expected.split("\n")
        if not isinstance(result, list):
            result = result.split("\n")
        for idx, (expected_line, result_line) in enumerate(zip(expected, result)):
            self.assertEqual(expected_line, result_line,
                             "Line %d:\nExp: %s\nGot: %s" % (idx, expected_line, result_line))
        self.assertEqual(len(expected), len(result),
                         "Unmatched lines. Got:\n%s\nExpected:\n%s" % ("\n".join(expected), "\n".join(result)))

    def codeToLines(self, tree):
        writer = CodeWriter()
        writer.write(tree)
        return writer.result.lines

    def codeToString(self, tree):
        return "\n".join(self.codeToLines(tree))

    def assertCode(self, expected, result_tree):
        result_lines = self.codeToLines(result_tree)

        expected_lines = strip_common_indent(expected.split("\n"))

        for idx, (line, expected_line) in enumerate(zip(result_lines, expected_lines)):
            self.assertEqual(expected_line, line,
                             "Line %d:\nGot: %s\nExp: %s" % (idx, line, expected_line))
        self.assertEqual(len(result_lines), len(expected_lines),
                         "Unmatched lines. Got:\n%s\nExpected:\n%s" % ("\n".join(result_lines), expected))

    def assertNodeExists(self, path, result_tree):
        self.assertNotEqual(TreePath.find_first(result_tree, path), None,
                            "Path '%s' not found in result tree" % path)

    def fragment(self, code, pxds=None, pipeline=None):
        "Simply create a tree fragment using the name of the test-case in parse errors."
        if pxds is None:
            pxds = {}
        if pipeline is None:
            pipeline = []
        name = self.id()
        if name.startswith("__main__."):
            name = name[len("__main__."):]
        name = name.replace(".", "_")
        return TreeFragment(code, name, pxds, pipeline=pipeline)

    def treetypes(self, root):
        return treetypes(root)

    def should_fail(self, func, exc_type=Exception):
        """Calls "func" and fails if it doesn't raise the right exception
        (any exception by default). Also returns the exception in question.
        """
        try:
            func()
            self.fail("Expected an exception of type %r" % exc_type)
        except exc_type as e:
            self.assertTrue(isinstance(e, exc_type))
            return e

    def should_not_fail(self, func):
        """Calls func and succeeds if and only if no exception is raised
        (i.e. converts exception raising into a failed testcase). Returns
        the return value of func."""
        try:
            return func()
        except Exception as exc:
            self.fail(str(exc))


class TransformTest(CythonTest):
    """
    Utility base class for transform unit tests. It is based around constructing
    test trees (either explicitly or by parsing a Cython code string); running
    the transform, serialize it using a customized Cython serializer (with
    special markup for nodes that cannot be represented in Cython),
    and do a string-comparison line-by-line of the result.

    To create a test case:
     - Call run_pipeline. The pipeline should at least contain the transform you
       are testing; pyx should be either a string (passed to the parser to
       create a post-parse tree) or a node representing input to pipeline.
       The result will be a transformed result.

     - Check that the tree is correct. If wanted, assertCode can be used, which
       takes a code string as expected, and a ModuleNode in result_tree
       (it serializes the ModuleNode to a string and compares line-by-line).

    All code strings are first stripped for whitespace lines and then common
    indentation.

    Plans: One could have a pxd dictionary parameter to run_pipeline.
    """

    def run_pipeline(self, pipeline, pyx, pxds=None):
        if pxds is None:
            pxds = {}
        tree = self.fragment(pyx, pxds).root
        # Run pipeline
        for T in pipeline:
            tree = T(tree)
        return tree


# For the test C code validation, we have to take care that the test directives (and thus
# the match strings) do not just appear in (multiline) C code comments containing the original
# Cython source code.  Thus, we discard the comments before matching.
# This seems a prime case for re.VERBOSE, but it seems to match some of the whitespace.
_strip_c_comments = partial(re.compile(
    re.sub(r'\s+', '', r'''
        /[*] (
            (?: [^*\n] | [*][^/] )*
            [\n]
            (?: [^*] | [*][^/] )*
        ) [*]/
    ''')
).sub, '')

_strip_cython_code_from_html = partial(re.compile(
    re.sub(r'\s\s+', '', r'''
    (?:
        <pre class=["'][^"']*cython\s+line[^"']*["']\s*>
        (?:[^<]|<(?!/pre))+
        </pre>
    )|(?:
        <style[^>]*>
        (?:[^<]|<(?!/style))+
        </style>
    )
    ''')
).sub, '')


def _parse_pattern(pattern):
    start = end = None
    if pattern.startswith('/'):
        start, pattern = re.split(r"(?<!\\)/", pattern[1:], maxsplit=1)
        pattern = pattern.strip()
    if pattern.startswith(':'):
        pattern = pattern[1:].strip()
        if pattern.startswith("/"):
            end, pattern = re.split(r"(?<!\\)/", pattern[1:], maxsplit=1)
            pattern = pattern.strip()
    return start, end, pattern


class TreeAssertVisitor(VisitorTransform):
    # actually, a TreeVisitor would be enough, but this needs to run
    # as part of the compiler pipeline

    def __init__(self):
        super().__init__()
        self._module_pos = None
        self._c_patterns = []
        self._c_antipatterns = []

    def create_c_file_validator(self):
        patterns, antipatterns = self._c_patterns, self._c_antipatterns

        def fail(pos, pattern, found, file_path):
            Errors.error(pos, "Pattern '%s' %s found in %s" %(
                pattern,
                'was' if found else 'was not',
                file_path,
            ))

        def extract_section(file_path, content, start, end):
            if start:
                split = re.search(start, content)
                if split:
                    content = content[split.end():]
                else:
                    fail(self._module_pos, start, found=False, file_path=file_path)
            if end:
                split = re.search(end, content)
                if split:
                    content = content[:split.start()]
                else:
                    fail(self._module_pos, end, found=False, file_path=file_path)
            return content

        def validate_file_content(file_path, content):
            for pattern in patterns:
                #print("Searching pattern '%s'" % pattern)
                start, end, pattern = _parse_pattern(pattern)
                section = extract_section(file_path, content, start, end)
                if not re.search(pattern, section):
                    fail(self._module_pos, pattern, found=False, file_path=file_path)

            for antipattern in antipatterns:
                #print("Searching antipattern '%s'" % antipattern)
                start, end, antipattern = _parse_pattern(antipattern)
                section = extract_section(file_path, content, start, end)
                if re.search(antipattern, section):
                    fail(self._module_pos, antipattern, found=True, file_path=file_path)

        def validate_c_file(result):
            c_file = result.c_file
            if not (patterns or antipatterns):
                #print("No patterns defined for %s" % c_file)
                return result

            with open(c_file, encoding='utf8') as f:
                content = f.read()
            content = _strip_c_comments(content)
            validate_file_content(c_file, content)

        return validate_c_file

    def _check_directives(self, node):
        directives = node.directives
        if 'test_assert_path_exists' in directives:
            for path in directives['test_assert_path_exists']:
                if TreePath.find_first(node, path) is None:
                    Errors.error(
                        node.pos,
                        "Expected path '%s' not found in result tree" % path)
        if 'test_fail_if_path_exists' in directives:
            for path in directives['test_fail_if_path_exists']:
                first_node = TreePath.find_first(node, path)
                if first_node is not None:
                    Errors.error(
                        first_node.pos,
                        "Unexpected path '%s' found in result tree" % path)
        if 'test_assert_c_code_has' in directives:
            self._c_patterns.extend(directives['test_assert_c_code_has'])
        if 'test_fail_if_c_code_has' in directives:
            self._c_antipatterns.extend(directives['test_fail_if_c_code_has'])

    def visit_ModuleNode(self, node):
        self._module_pos = node.pos
        self._check_directives(node)
        self.visitchildren(node)
        return node

    def visit_CompilerDirectivesNode(self, node):
        self._check_directives(node)
        self.visitchildren(node)
        return node

    visit_Node = VisitorTransform.recurse_to_children


def unpack_source_tree(tree_file, workdir, cython_root):
    programs = {
        'PYTHON': [sys.executable],
        'CYTHON': [sys.executable, os.path.join(cython_root, 'cython.py')],
        'CYTHONIZE': [sys.executable, os.path.join(cython_root, 'cythonize.py')]
    }

    if workdir is None:
        workdir = tempfile.mkdtemp()
    header, cur_file = [], None
    with open(tree_file, 'rb') as f:
        try:
            for line in f:
                if line[:5] == b'#####':
                    filename = line.strip().strip(b'#').strip().decode('utf8').replace('/', os.path.sep)
                    path = os.path.join(workdir, filename)
                    if not os.path.exists(os.path.dirname(path)):
                        os.makedirs(os.path.dirname(path))
                    if cur_file is not None:
                        to_close, cur_file = cur_file, None
                        to_close.close()
                    cur_file = open(path, 'wb')
                elif cur_file is not None:
                    cur_file.write(line)
                elif line.strip() and not line.lstrip().startswith(b'#'):
                    if line.strip() not in (b'"""', b"'''"):
                        command = shlex.split(line.decode('utf8'))
                        if not command: continue
                        # In Python 3: prog, *args = command
                        prog, args = command[0], command[1:]
                        try:
                            header.append(programs[prog]+args)
                        except KeyError:
                            header.append(command)
        finally:
            if cur_file is not None:
                cur_file.close()
    return workdir, header


def write_file(file_path, content, dedent=False, encoding=None):
    r"""Write some content (text or bytes) to the file
    at `file_path` without translating `'\n'` into `os.linesep`.

    The default encoding is `'utf-8'`.
    """
    if isinstance(content, bytes):
        mode = "wb"

        # binary mode doesn't take an encoding and newline arguments
        newline = None
        default_encoding = None
    else:
        mode = "w"

        # any "\n" characters written are not translated
        # to the system default line separator, os.linesep
        newline = "\n"
        default_encoding = "utf-8"

    if encoding is None:
        encoding = default_encoding

    if dedent:
        content = textwrap.dedent(content)

    with open(file_path, mode=mode, encoding=encoding, newline=newline) as f:
        f.write(content)


def write_newer_file(file_path, newer_than, content, dedent=False, encoding=None):
    r"""
    Write `content` to the file `file_path` without translating `'\n'`
    into `os.linesep` and make sure it is newer than the file `newer_than`.

    The default encoding is `'utf-8'` (same as for `write_file`).
    """
    write_file(file_path, content, dedent=dedent, encoding=encoding)

    try:
        other_time = os.path.getmtime(newer_than)
    except OSError:
        # Support writing a fresh file (which is always newer than a non-existent one)
        other_time = None

    while other_time is None or other_time >= os.path.getmtime(file_path):
        write_file(file_path, content, dedent=dedent, encoding=encoding)


def py_parse_code(code):
    """
    Compiles code far enough to get errors from the parser and post-parse stage.

    Is useful for checking for syntax errors, however it doesn't generate runable
    code.
    """
    context = StringParseContext("test")
    # all the errors we care about are in the parsing or postparse stage
    try:
        with Errors.local_errors() as errors:
            result = TreeFragment(code, pipeline=[PostParse(context)])
            result = result.substitute()
        if errors:
            raise errors[0]  # compile error, which should get caught below
        else:
            return result
    except Errors.CompileError as e:
        raise SyntaxError(e.message_only)
