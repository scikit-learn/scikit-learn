r"""Parser for python source files."""

# Created Sun Nov 27 14:03:07 2016
# Author: Óscar Nájera

import ast
import re
import tokenize
from collections import namedtuple
from io import BytesIO
from textwrap import dedent

from sphinx.errors import ExtensionError
from sphinx.util.logging import getLogger

logger = getLogger("sphinx-gallery")

SYNTAX_ERROR_DOCSTRING = """
SyntaxError
===========

Example script with invalid Python syntax
"""

# The pattern for in-file config comments is designed to not greedily match
# newlines at the start and end, except for one newline at the end. This
# ensures that the matched pattern can be removed from the code without
# changing the block structure; i.e. empty newlines are preserved, e.g. in
#
#     a = 1
#
#     # sphinx_gallery_thumbnail_number = 2
#
#     b = 2
FLAG_START = r"^[\ \t]*#\s*"
FLAG_BODY = r"sphinx_gallery_([A-Za-z0-9_]+)(\s*=\s*(.+))?[\ \t]*\n?"
INFILE_CONFIG_PATTERN = re.compile(FLAG_START + FLAG_BODY, re.MULTILINE)

START_IGNORE_FLAG = FLAG_START + "sphinx_gallery_start_ignore"
END_IGNORE_FLAG = FLAG_START + "sphinx_gallery_end_ignore"
IGNORE_BLOCK_PATTERN = re.compile(
    rf"{START_IGNORE_FLAG}(?:[\s\S]*?){END_IGNORE_FLAG}\n?", re.MULTILINE
)


def parse_source_file(filename):
    """Parse source file into AST node.

    Parameters
    ----------
    filename : str
        File path

    Returns
    -------
    node : AST node
    content : utf-8 encoded string
    """
    # builtin open automatically converts \r\n to \n
    with open(filename, "r", encoding="utf-8") as fid:
        content = fid.read()

    try:
        node = ast.parse(content)
        return node, content
    except SyntaxError:
        return None, content


def _get_docstring_and_rest(filename):
    """Separate ``filename`` content between docstring and the rest.

    Strongly inspired from ast.get_docstring.

    Returns
    -------
    docstring : str
        docstring of ``filename``
    rest : str
        ``filename`` content without the docstring
    lineno : int
        The line number.
    node : ast.Module
        The ast node. When `filename` parsed with `mode='exec'` node should be
        of type `ast.Module`.
    """
    node, content = parse_source_file(filename)

    if node is None:
        return SYNTAX_ERROR_DOCSTRING, content, 1, node

    if not isinstance(node, ast.Module):
        raise ExtensionError(
            "This function only supports modules. You provided {}".format(
                node.__class__.__name__
            )
        )
    if not (
        node.body
        and isinstance(node.body[0], ast.Expr)
        and isinstance(node.body[0].value, ast.Constant)
    ):
        raise ExtensionError(
            'Could not find docstring in file "{}". '
            "A docstring is required by sphinx-gallery "
            'unless the file is ignored by "ignore_pattern"'.format(filename)
        )

    # Python 3.7+ way
    docstring = ast.get_docstring(node)
    assert docstring is not None  # should be guaranteed above
    # This is just for backward compat
    if node.body[0].value.value[:1] == "\n":
        # just for strict backward compat here
        docstring = "\n" + docstring
    ts = tokenize.tokenize(BytesIO(content.encode()).readline)
    # find the first string according to the tokenizer and get its end row
    for tk in ts:
        if tk.exact_type == 3:
            lineno, _ = tk.end
            break
    else:
        lineno = 0

    # This get the content of the file after the docstring last line
    # Note: 'maxsplit' argument is not a keyword argument in python2
    rest = "\n".join(content.split("\n")[lineno:])
    lineno += 1
    return docstring, rest, lineno, node


def extract_file_config(content):
    """Pull out the file-specific config specified in the docstring."""
    file_conf = {}
    for match in re.finditer(INFILE_CONFIG_PATTERN, content):
        name = match.group(1)
        value = match.group(3)
        if value is None:  # a flag rather than a config setting
            continue
        try:
            value = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            logger.warning(
                "Sphinx-gallery option %s was passed invalid value %s", name, value
            )
        else:
            file_conf[name] = value
    return file_conf


Block = namedtuple("Block", ["type", "content", "lineno"])
# type: "text" or "code"
# content (str): the block lines as str
# lineno (int): the line number where the block starts


def split_code_and_text_blocks(source_file, return_node=False):
    """Return list with source file separated into code and text blocks.

    Parameters
    ----------
    source_file : str
        Path to the source file.
    return_node : bool
        If True, return the ast node.

    Returns
    -------
    file_conf : dict
        File-specific settings given in source file comments as:
        ``# sphinx_gallery_<name> = <value>``
    blocks : list
        (label, content, line_number)
        List where each element is a tuple with the label ('text' or 'code'),
        the corresponding content string of block and the leading line number.
    node : ast.Module
        The parsed ast node.
    """
    docstring, rest_of_content, lineno, node = _get_docstring_and_rest(source_file)
    blocks = [Block("text", docstring, 1)]

    file_conf = extract_file_config(rest_of_content)

    pattern = re.compile(
        r"(?P<header_line>^#{20,}.*|^# ?%%.*)\s(?P<text_content>(?:^#.*\s?)*)",
        flags=re.M,
    )
    sub_pat = re.compile("^#", flags=re.M)

    pos_so_far = 0
    for match in re.finditer(pattern, rest_of_content):
        code_block_content = rest_of_content[pos_so_far : match.start()]
        if code_block_content.strip():
            blocks.append(Block("code", code_block_content, lineno))
        lineno += code_block_content.count("\n")

        lineno += 1  # Ignored header line of hashes.
        text_content = match.group("text_content")
        text_block_content = dedent(re.sub(sub_pat, "", text_content)).lstrip()
        if text_block_content.strip():
            blocks.append(Block("text", text_block_content, lineno))
        lineno += text_content.count("\n")

        pos_so_far = match.end()

    remaining_content = rest_of_content[pos_so_far:]
    if remaining_content.strip():
        blocks.append(Block("code", remaining_content, lineno))

    out = (file_conf, blocks)
    if return_node:
        out += (node,)
    return out


def remove_ignore_blocks(code_block):
    """
    Return the content of *code_block* with ignored areas removed.

    An ignore block starts with # sphinx_gallery_start_ignore, and ends with
    # sphinx_gallery_end_ignore. These lines and anything in between them will
    be removed, but surrounding empty lines are preserved.

    Parameters
    ----------
    code_block : str
        A code segment.
    """
    num_start_flags = len(re.findall(START_IGNORE_FLAG, code_block, re.MULTILINE))
    num_end_flags = len(re.findall(END_IGNORE_FLAG, code_block, re.MULTILINE))

    if num_start_flags != num_end_flags:
        raise ExtensionError(
            'All "sphinx_gallery_start_ignore" flags must have a matching '
            '"sphinx_gallery_end_ignore" flag!'
        )
    return re.subn(IGNORE_BLOCK_PATTERN, "", code_block)[0]


def remove_config_comments(code_block):
    """
    Return the content of *code_block* with in-file config comments removed.

    Comment lines of the pattern '# sphinx_gallery_[option] = [val]' are
    removed, but surrounding empty lines are preserved.

    Parameters
    ----------
    code_block : str
        A code segment.
    """
    parsed_code, _ = re.subn(INFILE_CONFIG_PATTERN, "", code_block)
    return parsed_code
