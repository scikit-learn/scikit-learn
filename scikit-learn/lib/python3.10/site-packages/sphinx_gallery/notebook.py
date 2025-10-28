r"""Parser for Jupyter notebooks.

Class that holds the Jupyter notebook information
"""

# Author: Óscar Nájera
# License: 3-clause BSD

import argparse
import base64
import copy
import json
import mimetypes
import os
import re
import sys
import textwrap
from collections import defaultdict
from functools import partial
from itertools import count
from pathlib import Path

import sphinx.util
from sphinx.errors import ExtensionError

from .py_source_parser import split_code_and_text_blocks

logger = sphinx.util.logging.getLogger("sphinx-gallery")


def jupyter_notebook_skeleton():
    """Returns a dictionary with the elements of a Jupyter notebook."""
    py_version = sys.version_info
    notebook_skeleton = {
        "cells": [],
        "metadata": {
            "kernelspec": {
                "display_name": "Python " + str(py_version[0]),
                "language": "python",
                "name": "python" + str(py_version[0]),
            },
            "language_info": {
                "codemirror_mode": {"name": "ipython", "version": py_version[0]},
                "file_extension": ".py",
                "mimetype": "text/x-python",
                "name": "python",
                "nbconvert_exporter": "python",
                "pygments_lexer": "ipython" + str(py_version[0]),
                "version": "{}.{}.{}".format(*sys.version_info[:3]),
            },
        },
        "nbformat": 4,
        "nbformat_minor": 0,
    }
    return notebook_skeleton


def directive_fun(match, directive):
    """Helper to fill in directives."""
    directive_to_alert = dict(note="info", warning="danger")
    return '<div class="alert alert-{}"><h4>{}</h4><p>{}</p></div>'.format(
        directive_to_alert[directive], directive.capitalize(), match.group(1).strip()
    )


def convert_code_to_md(text):
    """Rewrite code blocks to use Markdown's preferred backtick notation.

    Backtick notation preserves syntax highlighting.

    Parameters
    ----------
    text: str
        A mostly converted string of markdown text. May contain zero, one,
        or multiple code blocks in code-block format.
    """
    code_regex = r"[ \t]*\.\. code-block::[ \t]*(\S*)\n[ \t]*\n([ \t]+)"
    while True:
        code_block = re.search(code_regex, text)
        if not code_block:
            break
        indent = code_block.group(2)
        start_index = code_block.span()[1] - len(indent)

        # Find first non-empty, non-indented line
        end = re.compile(rf"^(?!{re.escape(indent)})[ \t]*\S+", re.MULTILINE)
        code_end_match = end.search(text, start_index)
        end_index = code_end_match.start() if code_end_match else len(text)

        contents = textwrap.dedent(text[start_index:end_index]).rstrip()
        new_code = f"```{code_block.group(1)}\n{contents}\n```\n"
        text = text[: code_block.span()[0]] + new_code + text[end_index:]
    return text


def rst2md(text, gallery_conf, target_dir, heading_levels):
    """Converts reST from docstrings and text blocks to markdown for Jupyter notebooks.

    Parameters
    ----------
    text: str
        reST input to be converted to MD
    gallery_conf : dict
        The sphinx-gallery configuration dictionary.
    target_dir : str
        Path that notebook is intended for. Used where relative paths
        may be required.
    heading_levels: dict
        Mapping of heading style ``(over_char, under_char)`` to heading level.
        Note that ``over_char`` is `None` when only underline is present.
    """
    # Characters recommended for use with headings
    # https://docutils.readthedocs.io/en/sphinx-docs/user/rst/quickstart.html#sections
    adornment_characters = "=`:.'\"~^_*+#<>-"
    headings = re.compile(
        # Start of string or blank line
        r"(?P<pre>\A|^[ \t]*\n)"
        # Optional over characters, allowing leading space on heading text
        r"(?:(?P<over>[{0}])(?P=over)*\n[ \t]*)?"
        # The heading itself, with at least one non-white space character
        r"(?P<heading>\S[^\n]*)\n"
        # Under character, setting to same character if over present.
        r"(?P<under>(?(over)(?P=over)|[{0}]))(?P=under)*$"
        r"".format(adornment_characters),
        flags=re.M,
    )

    text = re.sub(
        headings,
        lambda match: "{1}{0} {2}".format(
            "#" * heading_levels[match.group("over", "under")],
            *match.group("pre", "heading"),
        ),
        text,
    )

    math_eq = re.compile(r"^\.\. math::((?:.+)?(?:\n+^  .+)*)", flags=re.M)
    text = re.sub(
        math_eq,
        lambda match: r"\begin{{align}}{0}\end{{align}}".format(match.group(1).strip()),
        text,
    )
    inline_math = re.compile(r":math:`(.+?)`", re.DOTALL)
    text = re.sub(inline_math, r"$\1$", text)

    directives = ("warning", "note")
    for directive in directives:
        directive_re = re.compile(
            r"^\.\. %s::((?:.+)?(?:\n+^  .+)*)" % directive, flags=re.M
        )
        text = re.sub(directive_re, partial(directive_fun, directive=directive), text)

    footnote_links = re.compile(r"^ *\.\. _.*:.*$\n", flags=re.M)
    text = re.sub(footnote_links, "", text)

    embedded_uris = re.compile(r"`([^`]*?)\s*<([^`]*)>`_")
    text = re.sub(embedded_uris, r"[\1](\2)", text)

    refs = re.compile(r":ref:`")
    text = re.sub(refs, "`", text)

    contents = re.compile(r"^\s*\.\. contents::.*$(\n +:\S+: *$)*\n", flags=re.M)
    text = re.sub(contents, "", text)

    images = re.compile(r"^\.\. image::(.*$)((?:\n +:\S+:.*$)*)\n", flags=re.M)
    image_opts = re.compile(r"\n +:(\S+): +(.*)$", flags=re.M)
    text = re.sub(
        images,
        lambda match: '<img src="{}"{}>\n'.format(
            generate_image_src(match.group(1).strip(), gallery_conf, target_dir),
            re.sub(image_opts, r' \1="\2"', match.group(2) or ""),
        ),
        text,
    )

    text = convert_code_to_md(text)

    return text


def generate_image_src(image_path, gallery_conf, target_dir):
    """Modify image path for notebook, according to "notebook_images" config.

    URLs are unchanged.
    If "notebook_images" config is a str, it is used as a prefix to image path, relative
    to "src_dir". If "notebook_images" is `True`, image is embedded as URI. If
    "notebook_images" is `False`, "file://" is prepended.
    """
    if re.match(r"https?://", image_path):
        return image_path

    if not gallery_conf["notebook_images"]:
        return "file://" + image_path.lstrip("/")

    # If absolute path from source directory given
    if image_path.startswith("/"):
        # Path should now be relative to source dir, not target dir
        target_dir = gallery_conf["src_dir"]
        image_path = image_path.lstrip("/")
    full_path = os.path.join(target_dir, image_path.replace("/", os.sep))

    if isinstance(gallery_conf["notebook_images"], str):
        # Use as prefix e.g. URL
        prefix = gallery_conf["notebook_images"]
        rel_path = os.path.relpath(full_path, gallery_conf["src_dir"])
        return prefix + rel_path.replace(os.sep, "/")
    else:
        # True, but not string. Embed as data URI.
        try:
            with open(full_path, "rb") as image_file:
                data = base64.b64encode(image_file.read())
        except OSError:
            raise ExtensionError(
                f"Unable to open {full_path} to generate notebook data URI"
            )
        mime_type = mimetypes.guess_type(full_path)
        return f"data:{mime_type[0]};base64,{data.decode('ascii')}"


def jupyter_notebook(script_blocks, gallery_conf, target_dir):
    """Generate a Jupyter notebook file cell-by-cell.

    Parameters
    ----------
    script_blocks : list
        Script execution cells.
    gallery_conf : dict
        The sphinx-gallery configuration dictionary.
    target_dir : str
        Path that notebook is intended for. Used where relative paths
        may be required.
    """
    first_cell = gallery_conf["first_notebook_cell"]
    last_cell = gallery_conf["last_notebook_cell"]
    work_notebook = jupyter_notebook_skeleton()
    if first_cell is not None:
        add_code_cell(work_notebook, first_cell)
    fill_notebook(work_notebook, script_blocks, gallery_conf, target_dir)
    if last_cell is not None:
        add_code_cell(work_notebook, last_cell)

    return work_notebook


def add_code_cell(work_notebook, code):
    """Add a code cell to the notebook.

    Parameters
    ----------
    code : str
        Cell content
    """
    code_cell = {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {"collapsed": False},
        "outputs": [],
        "source": [code.strip()],
    }
    work_notebook["cells"].append(code_cell)


def add_markdown_cell(work_notebook, markdown):
    """Add a markdown cell to the notebook.

    Parameters
    ----------
    markdown : str
        Markdown cell content.
    """
    markdown_cell = {"cell_type": "markdown", "metadata": {}, "source": [markdown]}
    work_notebook["cells"].append(markdown_cell)


def promote_jupyter_cell_magic(work_notebook, markdown):
    """Promote Jupyter cell magic in text blocks to code block in notebooks.

    Parses a block of markdown text looking for code blocks starting with a
    Jupyter cell magic (e.g. %%bash). Whenever one is found, the text before it
    and the code (as a runnable code block) are added to work_notebook. Any
    remaining text is returned.

    Parameters
    ----------
    markdown : str
        Markdown cell content.
    """
    # Regex detects all code blocks that use %% Jupyter cell magic
    cell_magic_regex = r"\n?```\s*[a-z]*\n(%%(?:[\s\S]*?))\n?```\n?"

    text_cell_start = 0
    for magic_cell in re.finditer(cell_magic_regex, markdown):
        # Extract the preceding text block, and add it if non-empty
        text_block = markdown[text_cell_start : magic_cell.span()[0]]
        if text_block and not text_block.isspace():
            add_markdown_cell(work_notebook, text_block)
        text_cell_start = magic_cell.span()[1]

        code_block = magic_cell.group(1)
        add_code_cell(work_notebook, code_block)

    # Return remaining text (which equals markdown if no magic cells exist)
    return markdown[text_cell_start:]


def fill_notebook(work_notebook, script_blocks, gallery_conf, target_dir):
    """Writes the Jupyter notebook cells.

    If available, uses pypandoc to convert rst to markdown.

    Parameters
    ----------
    script_blocks : list
        Each list element should be a tuple of (label, content, lineno).
    """
    heading_level_counter = count(start=1)
    heading_levels = defaultdict(lambda: next(heading_level_counter))
    for blabel, bcontent, lineno in script_blocks:
        if blabel == "code":
            add_code_cell(work_notebook, bcontent)
        else:
            if gallery_conf["pypandoc"] is False:
                markdown = rst2md(
                    bcontent + "\n", gallery_conf, target_dir, heading_levels
                )
            else:
                import pypandoc

                # pandoc automatically adds \n to the end
                markdown = pypandoc.convert_text(
                    bcontent, to="md", format="rst", **gallery_conf["pypandoc"]
                )

            if gallery_conf["promote_jupyter_magic"]:
                remaining = promote_jupyter_cell_magic(work_notebook, markdown)
                if remaining and not remaining.isspace():
                    add_markdown_cell(work_notebook, remaining)
            else:
                add_markdown_cell(work_notebook, markdown)


def save_notebook(work_notebook, write_file):
    """Saves the Jupyter work_notebook to write_file."""
    with open(write_file, "w") as out_nb:
        json.dump(work_notebook, out_nb, indent=2)


###############################################################################
# Notebook shell utility


def python_to_jupyter_cli(args=None, namespace=None, sphinx_gallery_conf=None):
    """Exposes the jupyter notebook renderer to the command line.

    Takes the same arguments as ArgumentParser.parse_args.
    `sphinx_gallery_conf` functions same as in `conf.py`.
    """
    from . import gen_gallery  # To avoid circular import

    parser = argparse.ArgumentParser(description="Sphinx-Gallery Notebook converter")
    parser.add_argument(
        "python_src_file",
        nargs="+",
        help="Input Python file script to convert. "
        "Supports multiple files and shell wildcards"
        " (e.g. *.py)",
    )
    args = parser.parse_args(args, namespace)

    # handle `sphinx_gallery_conf`
    gallery_conf = copy.deepcopy(gen_gallery.DEFAULT_GALLERY_CONF)
    if sphinx_gallery_conf is not None:
        gallery_conf.update(sphinx_gallery_conf)

    # run script
    for src_file in args.python_src_file:
        file_conf, blocks = split_code_and_text_blocks(src_file)
        print(f"Converting {src_file}")
        target_dir = os.path.dirname(src_file)
        example_nb = jupyter_notebook(blocks, copy.deepcopy(gallery_conf), target_dir)
        save_notebook(example_nb, Path(src_file).with_suffix(".ipynb"))
