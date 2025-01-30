"""BlockParser divides non `.py` source files into blocks of code and markup text."""

import ast
import codecs
import re
from pathlib import Path
from textwrap import dedent

import pygments.lexers
import pygments.token
from sphinx.errors import ExtensionError
from sphinx.util.logging import getLogger

from .py_source_parser import FLAG_BODY, Block

logger = getLogger("sphinx-gallery")

# Don't just use "x in pygments.token.Comment" because it also includes preprocessor
# statements
COMMENT_TYPES = (
    pygments.token.Comment.Single,
    pygments.token.Comment,
    pygments.token.Comment.Multiline,
)


class BlockParser:
    """
    A parser that breaks a source file into blocks of code and markup text.

    Determines the source language and identifies comment blocks using pygments.

    Parameters
    ----------
    source_file : str
        A file name that has a suffix compatible with files that are subsequently parsed
    gallery_conf : dict
        Contains the configuration of Sphinx-Gallery.
    """

    def __init__(self, source_file, gallery_conf):
        source_file = Path(source_file)
        if name := gallery_conf["filetype_parsers"].get(source_file.suffix):
            self.lexer = pygments.lexers.find_lexer_class_by_name(name)()
        else:
            self.lexer = pygments.lexers.find_lexer_class_for_filename(source_file)()
        self.language = self.lexer.name

        # determine valid comment syntaxes. For each possible syntax, the tuple contains
        # - A test comment
        # - The comment start character
        # - The comment end character (for multiline comments) or None
        # - A regex for the start of a block based on repeated characters, or None
        comment_tests = [
            ("#= comment =#", "#=", "=#", None),  # Julia multiline
            ("# comment", "#", None, "#{20,}"),  # Julia, Ruby, Bash, Perl, etc.
            ("// comment", "//", None, "/{20,}"),  # C++, C#, Java, Rust, etc.
            ("/* comment */", r"/\*", r"\*/", r"/\*{20,}/"),  # C/C++ etc. multiline
            ("% comment", "%", None, "%{20,}"),  # Matlab
            ("! comment", "!", None, "!{20,}"),  # Fortran
            ("c     comment", r"^c(?:$|     )", None, None),  # Fortran 77
        ]

        self.allowed_comments = []
        allowed_special = []
        self.multiline_end = re.compile(chr(0))  # unmatchable regex
        for test, start, end, special in comment_tests:
            if next(self.lexer.get_tokens(test))[0] in COMMENT_TYPES:
                self.allowed_comments.append(start)
                if end:
                    self.multiline_end = re.compile(rf"(.*?)\s*{end}")
                if special:
                    allowed_special.append(special)

        if self.language == "Matlab":
            # Matlab treats code sections starting with "%%" in a similar manner to
            # tools that recognize "# %%" in Python code, so we accept that style as
            # an alternative.
            allowed_special.append(r"%%(?:$|\s)")

        if r"/\*" in self.allowed_comments:
            # Remove decorative asterisks and comment starts from C-style multiline
            # comments
            self.multiline_cleanup = re.compile(r"\s*/?\*\s*")
        else:
            self.multiline_cleanup = re.compile(r"\s*")

        comment_start = "|".join(self.allowed_comments)
        allowed_special = "|".join(allowed_special)
        if allowed_special:
            self.start_special = re.compile(
                f"(?:(?:{comment_start}) ?%% ?|{allowed_special})(.*)"
            )
        else:
            self.start_special = re.compile(f"(?:{comment_start}) ?%% ?(.*)")
        self.continue_text = re.compile(f"(?:{comment_start}) ?(.*)")

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
        flag_start = rf"^[\ \t]*(?:{comment_start})\s*"

        self.infile_config_pattern = re.compile(flag_start + FLAG_BODY, re.MULTILINE)
        self.start_ignore_flag = flag_start + "sphinx_gallery_start_ignore"
        self.end_ignore_flag = flag_start + "sphinx_gallery_end_ignore"
        self.ignore_block_pattern = re.compile(
            rf"{self.start_ignore_flag}(?:[\s\S]*?){self.end_ignore_flag}\n?",
            re.MULTILINE,
        )

    def split_code_and_text_blocks(self, source_file, return_node=False):
        """Return list with source file separated into code and text blocks.

        Parameters
        ----------
        source_file : str
            Path to the source file.
        return_node : bool
            Ignored; returning an ast node is not supported

        Returns
        -------
        file_conf : dict
            File-specific settings given in source file comments as:
            ``# sphinx_gallery_<name> = <value>``
        blocks : list
            (label, content, line_number)
            List where each element is a tuple with the label ('text' or 'code'),
            the corresponding content string of block and the leading line number.
        node : None
            Returning an ast node is not supported.
        """
        with codecs.open(source_file, "r", "utf-8") as fid:
            content = fid.read()
        # change from Windows format to UNIX for uniformity
        content = content.replace("\r\n", "\n")
        return self._split_content(content)

    def _get_content_lines(self, content):
        """
        Combine individual tokens into lines.

        Use the first non-whitespace token (if any) as the characteristic token type for
        the line.
        """
        current_line = []
        line_token = pygments.token.Whitespace
        for token, text in self.lexer.get_tokens(content):
            if line_token == pygments.token.Whitespace:
                line_token = token

            if "\n" in text:
                text_lines = text.split("\n")
                # first item belongs to the previous line
                current_line.append(text_lines.pop(0))
                yield line_token, "".join(current_line)
                # last item belongs to the line after this token
                current_line = [text_lines.pop()]
                # Anything left is a block of lines to add directly
                for ln in text_lines:
                    if not ln.strip():
                        line_token = pygments.token.Whitespace
                    yield line_token, ln
                if not current_line[0].strip():
                    line_token = pygments.token.Whitespace
            else:
                current_line.append(text)

    def _get_blocks(self, content):
        """
        Generate a sequence of "blocks" from the lines in ``content``.

        Each block is a tuple of (label, content, line_number), matching the format
        ultimately returned by `split_code_and_text_blocks`.
        """
        start_text = self.continue_text  # No special delimiter needed for first block
        needs_multiline_cleanup = False

        def cleanup_multiline(lines):
            nonlocal needs_multiline_cleanup
            first_line = 1 if start_text == self.continue_text else 0
            longest = max(len(line) for line in lines)
            matched = False
            for i, line in enumerate(lines[first_line:]):
                if m := self.multiline_cleanup.match(line):
                    matched = True
                    if (n := len(m.group(0))) < len(line):
                        longest = min(longest, n)

            if matched and longest:
                for i, line in enumerate(lines[first_line:], start=first_line):
                    lines[i] = lines[i][longest:]
            needs_multiline_cleanup = False
            return lines

        def finalize_block(mode, block):
            nonlocal start_text
            if mode == "text":
                if needs_multiline_cleanup:
                    cleanup_multiline(block)

                # subsequent blocks need to have the special delimiter
                start_text = self.start_special

                # Remove leading blank lines, and end in a single newline
                first = 0
                for i, line in enumerate(block):
                    first = i
                    if line.strip():
                        break
                last = None
                for i, line in enumerate(reversed(block)):
                    last = -i or None
                    if line.strip():
                        break
                block.append("")
                text = dedent("\n".join(block[first:last]))
            else:
                text = "\n".join(block)
            return Block(mode, text, n - len(block))

        block = []
        mode = None
        for n, (token, text) in enumerate(self._get_content_lines(content)):
            if mode == "text" and token in pygments.token.Whitespace:
                # Blank line ends current text block
                if block:
                    yield finalize_block(mode, block)
                mode, block = None, []
            elif (
                mode != "text"
                and token in COMMENT_TYPES
                and (m := start_text.search(text))
            ):
                # start of a text block; end the current block
                if block:
                    yield finalize_block(mode, block)
                mode, block = "text", []
                if (trailing_text := m.group(1)) is not None:
                    if start_text == self.continue_text:
                        if (delimited := self.start_special.search(text)) and (
                            heading := delimited.group(1).strip()
                        ):
                            # Treat text on the same line as the delimiter as a title
                            block.extend((heading, "=" * len(heading), ""))
                        else:
                            # Keep any text on the first line of the title block as is
                            block.append(trailing_text)
                    elif heading := trailing_text.strip():
                        # Treat text on the same line as the delimiter as a heading
                        block.extend((heading, "-" * len(heading), ""))
            elif mode == "text" and token in COMMENT_TYPES:
                # Continuation of a text block
                if m := self.start_special.search(text):
                    # Starting a new text block now is just a continuation of the
                    # existing one, possibly including a new heading.
                    if heading := m.group(1).strip():
                        block.extend((heading, "-" * len(heading), ""))
                    pass
                elif token == pygments.token.Comment.Multiline:
                    if m := self.multiline_end.search(text):
                        block.append(m.group(1))
                        needs_multiline_cleanup = True
                    else:
                        block.append(text)
                else:
                    block.append(self.continue_text.search(text).group(1))
            elif mode != "code":
                # start of a code block
                if block:
                    yield finalize_block(mode, block)
                mode, block = "code", [text]
            else:
                # continuation of a code block
                block.append(text)

        # end of input ends final block
        if block:
            yield finalize_block(mode, block)

    def _split_content(self, content):
        """
        Split the input content into blocks.

        Return a tuple of (file_conf, blocks, None) that corresponds to the return
        values of ``split_code_and_text_blocks``.
        """
        file_conf = self.extract_file_config(content)
        blocks = list(self._get_blocks(content))

        # For examples that start with a code block before the file docstring due to
        # language conventions or requirements, swap these blocks so the title block
        # comes first and merge consecutive code blocks if needed.
        if len(blocks) >= 2 and blocks[0].type == "code" and blocks[1].type == "text":
            blocks[0], blocks[1] = blocks[1], blocks[0]
            if len(blocks) >= 3 and blocks[2].type == "code":
                blocks[1] = Block(
                    "code",
                    f"{blocks[1].content}\n{blocks[2].content}",
                    blocks[1].lineno,
                )
                blocks.pop(2)

        return file_conf, blocks, None

    def extract_file_config(self, content):
        """Pull out the file-specific config specified in the docstring."""
        file_conf = {}
        for match in re.finditer(self.infile_config_pattern, content):
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

    def remove_ignore_blocks(self, code_block):
        """
        Return the content of *code_block* with ignored areas removed.

        An ignore block starts with ``?? sphinx_gallery_start_ignore`` and ends with
        ``?? sphinx_gallery_end_ignore`` where ``??`` is the active language's line
        comment marker. These lines and anything in between them will be removed, but
        surrounding empty lines are preserved.

        Parameters
        ----------
        code_block : str
            A code segment.
        """
        num_start_flags = len(re.findall(self.start_ignore_flag, code_block))
        num_end_flags = len(re.findall(self.end_ignore_flag, code_block))

        if num_start_flags != num_end_flags:
            raise ExtensionError(
                'All "sphinx_gallery_start_ignore" flags must have a matching '
                '"sphinx_gallery_end_ignore" flag!'
            )
        return re.subn(self.ignore_block_pattern, "", code_block)[0]

    def remove_config_comments(self, code_block):
        """
        Return the content of *code_block* with in-file config comments removed.

        Comment lines with the pattern ``sphinx_gallery_[option] = [val]`` after the
        line comment character are removed, but surrounding empty lines are preserved.

        Parameters
        ----------
        code_block : str
            A code segment.
        """
        parsed_code, _ = re.subn(self.infile_config_pattern, "", code_block)
        return parsed_code
