from __future__ import division, absolute_import, print_function

import sys
if sys.version_info[0] >= 3:
    from io import StringIO
else:
    from io import StringIO

import compiler
import inspect
import textwrap
import tokenize

from .compiler_unparse import unparse


class Comment(object):
    """ A comment block.
    """
    is_comment = True
    def __init__(self, start_lineno, end_lineno, text):
        # int : The first line number in the block. 1-indexed.
        self.start_lineno = start_lineno
        # int : The last line number. Inclusive!
        self.end_lineno = end_lineno
        # str : The text block including '#' character but not any leading spaces.
        self.text = text

    def add(self, string, start, end, line):
        """ Add a new comment line.
        """
        self.start_lineno = min(self.start_lineno, start[0])
        self.end_lineno = max(self.end_lineno, end[0])
        self.text += string

    def __repr__(self):
        return '%s(%r, %r, %r)' % (self.__class__.__name__, self.start_lineno,
            self.end_lineno, self.text)


class NonComment(object):
    """ A non-comment block of code.
    """
    is_comment = False
    def __init__(self, start_lineno, end_lineno):
        self.start_lineno = start_lineno
        self.end_lineno = end_lineno

    def add(self, string, start, end, line):
        """ Add lines to the block.
        """
        if string.strip():
            # Only add if not entirely whitespace.
            self.start_lineno = min(self.start_lineno, start[0])
            self.end_lineno = max(self.end_lineno, end[0])

    def __repr__(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self.start_lineno,
            self.end_lineno)


class CommentBlocker(object):
    """ Pull out contiguous comment blocks.
    """
    def __init__(self):
        # Start with a dummy.
        self.current_block = NonComment(0, 0)

        # All of the blocks seen so far.
        self.blocks = []

        # The index mapping lines of code to their associated comment blocks.
        self.index = {}

    def process_file(self, file):
        """ Process a file object.
        """
        if sys.version_info[0] >= 3:
            nxt = file.__next__
        else:
            nxt = file.next
        for token in tokenize.generate_tokens(nxt):
            self.process_token(*token)
        self.make_index()

    def process_token(self, kind, string, start, end, line):
        """ Process a single token.
        """
        if self.current_block.is_comment:
            if kind == tokenize.COMMENT:
                self.current_block.add(string, start, end, line)
            else:
                self.new_noncomment(start[0], end[0])
        else:
            if kind == tokenize.COMMENT:
                self.new_comment(string, start, end, line)
            else:
                self.current_block.add(string, start, end, line)

    def new_noncomment(self, start_lineno, end_lineno):
        """ We are transitioning from a noncomment to a comment.
        """
        block = NonComment(start_lineno, end_lineno)
        self.blocks.append(block)
        self.current_block = block

    def new_comment(self, string, start, end, line):
        """ Possibly add a new comment.

        Only adds a new comment if this comment is the only thing on the line.
        Otherwise, it extends the noncomment block.
        """
        prefix = line[:start[1]]
        if prefix.strip():
            # Oops! Trailing comment, not a comment block.
            self.current_block.add(string, start, end, line)
        else:
            # A comment block.
            block = Comment(start[0], end[0], string)
            self.blocks.append(block)
            self.current_block = block

    def make_index(self):
        """ Make the index mapping lines of actual code to their associated
        prefix comments.
        """
        for prev, block in zip(self.blocks[:-1], self.blocks[1:]):
            if not block.is_comment:
                self.index[block.start_lineno] = prev

    def search_for_comment(self, lineno, default=None):
        """ Find the comment block just before the given line number.

        Returns None (or the specified default) if there is no such block.
        """
        if not self.index:
            self.make_index()
        block = self.index.get(lineno, None)
        text = getattr(block, 'text', default)
        return text


def strip_comment_marker(text):
    """ Strip # markers at the front of a block of comment text.
    """
    lines = []
    for line in text.splitlines():
        lines.append(line.lstrip('#'))
    text = textwrap.dedent('\n'.join(lines))
    return text


def get_class_traits(klass):
    """ Yield all of the documentation for trait definitions on a class object.
    """
    # FIXME: gracefully handle errors here or in the caller?
    source = inspect.getsource(klass)
    cb = CommentBlocker()
    cb.process_file(StringIO(source))
    mod_ast = compiler.parse(source)
    class_ast = mod_ast.node.nodes[0]
    for node in class_ast.code.nodes:
        # FIXME: handle other kinds of assignments?
        if isinstance(node, compiler.ast.Assign):
            name = node.nodes[0].name
            rhs = unparse(node.expr).strip()
            doc = strip_comment_marker(cb.search_for_comment(node.lineno, default=''))
            yield name, rhs, doc

