# Copyright 2004-2005 Elemental Security, Inc. All Rights Reserved.
# Licensed to PSF under a Contributor Agreement.

# Modifications:
# Copyright 2006 Google, Inc. All Rights Reserved.
# Licensed to PSF under a Contributor Agreement.

"""Parser driver.

This provides a high-level interface to parse a file into a syntax tree.

"""

__author__ = "Guido van Rossum <guido@python.org>"

__all__ = ["Driver", "load_grammar"]

# Python imports
import io
import os
import logging
import pkgutil
import sys
from typing import (
    Any,
    IO,
    Iterable,
    List,
    Optional,
    Text,
    Tuple,
    Union,
)

# Pgen imports
from . import grammar, parse, token, tokenize, pgen
from logging import Logger
from blib2to3.pytree import _Convert, NL
from blib2to3.pgen2.grammar import Grammar

Path = Union[str, "os.PathLike[str]"]


class Driver(object):
    def __init__(
        self,
        grammar: Grammar,
        convert: Optional[_Convert] = None,
        logger: Optional[Logger] = None,
    ) -> None:
        self.grammar = grammar
        if logger is None:
            logger = logging.getLogger(__name__)
        self.logger = logger
        self.convert = convert

    def parse_tokens(self, tokens: Iterable[Any], debug: bool = False) -> NL:
        """Parse a series of tokens and return the syntax tree."""
        # XXX Move the prefix computation into a wrapper around tokenize.
        p = parse.Parser(self.grammar, self.convert)
        p.setup()
        lineno = 1
        column = 0
        indent_columns = []
        type = value = start = end = line_text = None
        prefix = ""
        for quintuple in tokens:
            type, value, start, end, line_text = quintuple
            if start != (lineno, column):
                assert (lineno, column) <= start, ((lineno, column), start)
                s_lineno, s_column = start
                if lineno < s_lineno:
                    prefix += "\n" * (s_lineno - lineno)
                    lineno = s_lineno
                    column = 0
                if column < s_column:
                    prefix += line_text[column:s_column]
                    column = s_column
            if type in (tokenize.COMMENT, tokenize.NL):
                prefix += value
                lineno, column = end
                if value.endswith("\n"):
                    lineno += 1
                    column = 0
                continue
            if type == token.OP:
                type = grammar.opmap[value]
            if debug:
                self.logger.debug(
                    "%s %r (prefix=%r)", token.tok_name[type], value, prefix
                )
            if type == token.INDENT:
                indent_columns.append(len(value))
                _prefix = prefix + value
                prefix = ""
                value = ""
            elif type == token.DEDENT:
                _indent_col = indent_columns.pop()
                prefix, _prefix = self._partially_consume_prefix(prefix, _indent_col)
            if p.addtoken(type, value, (prefix, start)):
                if debug:
                    self.logger.debug("Stop.")
                break
            prefix = ""
            if type in {token.INDENT, token.DEDENT}:
                prefix = _prefix
            lineno, column = end
            if value.endswith("\n"):
                lineno += 1
                column = 0
        else:
            # We never broke out -- EOF is too soon (how can this happen???)
            assert start is not None
            raise parse.ParseError("incomplete input", type, value, (prefix, start))
        assert p.rootnode is not None
        return p.rootnode

    def parse_stream_raw(self, stream: IO[Text], debug: bool = False) -> NL:
        """Parse a stream and return the syntax tree."""
        tokens = tokenize.generate_tokens(stream.readline, grammar=self.grammar)
        return self.parse_tokens(tokens, debug)

    def parse_stream(self, stream: IO[Text], debug: bool = False) -> NL:
        """Parse a stream and return the syntax tree."""
        return self.parse_stream_raw(stream, debug)

    def parse_file(
        self, filename: Path, encoding: Optional[Text] = None, debug: bool = False
    ) -> NL:
        """Parse a file and return the syntax tree."""
        with io.open(filename, "r", encoding=encoding) as stream:
            return self.parse_stream(stream, debug)

    def parse_string(self, text: Text, debug: bool = False) -> NL:
        """Parse a string and return the syntax tree."""
        tokens = tokenize.generate_tokens(
            io.StringIO(text).readline, grammar=self.grammar
        )
        return self.parse_tokens(tokens, debug)

    def _partially_consume_prefix(self, prefix: Text, column: int) -> Tuple[Text, Text]:
        lines: List[str] = []
        current_line = ""
        current_column = 0
        wait_for_nl = False
        for char in prefix:
            current_line += char
            if wait_for_nl:
                if char == "\n":
                    if current_line.strip() and current_column < column:
                        res = "".join(lines)
                        return res, prefix[len(res) :]

                    lines.append(current_line)
                    current_line = ""
                    current_column = 0
                    wait_for_nl = False
            elif char in " \t":
                current_column += 1
            elif char == "\n":
                # unexpected empty line
                current_column = 0
            else:
                # indent is finished
                wait_for_nl = True
        return "".join(lines), current_line


def _generate_pickle_name(gt: Path, cache_dir: Optional[Path] = None) -> Text:
    head, tail = os.path.splitext(gt)
    if tail == ".txt":
        tail = ""
    name = head + tail + ".".join(map(str, sys.version_info)) + ".pickle"
    if cache_dir:
        return os.path.join(cache_dir, os.path.basename(name))
    else:
        return name


def load_grammar(
    gt: Text = "Grammar.txt",
    gp: Optional[Text] = None,
    save: bool = True,
    force: bool = False,
    logger: Optional[Logger] = None,
) -> Grammar:
    """Load the grammar (maybe from a pickle)."""
    if logger is None:
        logger = logging.getLogger(__name__)
    gp = _generate_pickle_name(gt) if gp is None else gp
    if force or not _newer(gp, gt):
        logger.info("Generating grammar tables from %s", gt)
        g: grammar.Grammar = pgen.generate_grammar(gt)
        if save:
            logger.info("Writing grammar tables to %s", gp)
            try:
                g.dump(gp)
            except OSError as e:
                logger.info("Writing failed: %s", e)
    else:
        g = grammar.Grammar()
        g.load(gp)
    return g


def _newer(a: Text, b: Text) -> bool:
    """Inquire whether file a was written since file b."""
    if not os.path.exists(a):
        return False
    if not os.path.exists(b):
        return True
    return os.path.getmtime(a) >= os.path.getmtime(b)


def load_packaged_grammar(
    package: str, grammar_source: Text, cache_dir: Optional[Path] = None
) -> grammar.Grammar:
    """Normally, loads a pickled grammar by doing
        pkgutil.get_data(package, pickled_grammar)
    where *pickled_grammar* is computed from *grammar_source* by adding the
    Python version and using a ``.pickle`` extension.

    However, if *grammar_source* is an extant file, load_grammar(grammar_source)
    is called instead. This facilitates using a packaged grammar file when needed
    but preserves load_grammar's automatic regeneration behavior when possible.

    """
    if os.path.isfile(grammar_source):
        gp = _generate_pickle_name(grammar_source, cache_dir) if cache_dir else None
        return load_grammar(grammar_source, gp=gp)
    pickled_name = _generate_pickle_name(os.path.basename(grammar_source), cache_dir)
    data = pkgutil.get_data(package, pickled_name)
    assert data is not None
    g = grammar.Grammar()
    g.loads(data)
    return g


def main(*args: Text) -> bool:
    """Main program, when run as a script: produce grammar pickle files.

    Calls load_grammar for each argument, a path to a grammar text file.
    """
    if not args:
        args = tuple(sys.argv[1:])
    logging.basicConfig(level=logging.INFO, stream=sys.stdout, format="%(message)s")
    for gt in args:
        load_grammar(gt, save=True, force=True)
    return True


if __name__ == "__main__":
    sys.exit(int(not main()))
