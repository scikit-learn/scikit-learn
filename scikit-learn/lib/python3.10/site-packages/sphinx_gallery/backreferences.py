# Author: Óscar Nájera
# License: 3-clause BSD
"""Backreferences Generator.

Parses example file code in order to keep track of used functions.
"""

import ast
import inspect
import os
import re
import sys
from collections import defaultdict
from html import escape
from pathlib import Path

import sphinx.util
from sphinx.errors import ExtensionError

from ._dummy import DummyClass  # noqa: F401
from .scrapers import _find_image_ext
from .utils import _W_KW, _replace_md5

THUMBNAIL_PARENT_DIV = """
.. raw:: html

    <div class="sphx-glr-thumbnails">

.. thumbnail-parent-div-open
"""

THUMBNAIL_PARENT_DIV_CLOSE = """
.. thumbnail-parent-div-close

.. raw:: html

    </div>

"""


class NameFinder(ast.NodeVisitor):
    """Finds the longest form of variable names and their imports in code.

    Only retains names from imported modules.
    """

    def __init__(self, global_variables=None):
        super().__init__()
        self.imported_names = {}
        self.global_variables = global_variables or {}
        self.accessed_names = set()

    def visit_Import(self, node, prefix=""):
        """For 'import' add node names to `imported_names`."""
        for alias in node.names:
            local_name = alias.asname or alias.name
            self.imported_names[local_name] = prefix + alias.name

    def visit_ImportFrom(self, node):
        """For 'from import' add node names to `imported_names`, incl module prefix."""
        self.visit_Import(node, node.module + ".")

    def visit_Name(self, node):
        """Add node id to `accessed_names`."""
        self.accessed_names.add(node.id)

    def visit_Attribute(self, node):
        """Add attributes, including their prefix, to `accessed_names`."""
        attrs = []
        while isinstance(node, ast.Attribute):
            attrs.append(node.attr)
            node = node.value

        if isinstance(node, ast.Name):
            # This is a.b, not e.g. a().b
            attrs.append(node.id)
            self.accessed_names.add(".".join(reversed(attrs)))
        else:
            # need to get a in a().b
            self.visit(node)

    def get_mapping(self):
        """Map names used in code, using AST nodes, to their fully qualified names.

        Returns
        -------
        options : List[Tuple[str]]
            List of tuples, each tuple containing the following information about
            an `accessed_name`:

                accessed name : str,
                fully qualified name : str
                if it is a class attribute (i.e., property or method) : bool
                if it is a class : bool
                if it is an explicit backreference : bool (always false here)
        """
        options = list()
        for name in self.accessed_names:
            local_name_split = name.split(".")
            # first pass: by global variables and object inspection (preferred)
            for split_level in range(len(local_name_split)):
                local_name = ".".join(local_name_split[: split_level + 1])
                remainder = name[len(local_name) :]
                if local_name in self.global_variables:
                    obj = self.global_variables[local_name]
                    class_attr, method = False, []
                    if remainder:
                        for level in remainder[1:].split("."):
                            last_obj = obj
                            # determine if it's a property
                            prop = getattr(last_obj.__class__, level, None)
                            if isinstance(prop, property):
                                obj = last_obj
                                class_attr, method = True, [level]
                                break
                            # For most objects this will emit a AttributeError,
                            # but for some (e.g., PyQt5) you can get other
                            # errors like "RuntimeError: wrapped C/C++ object
                            # ... has been deleted" so let's be safer with
                            # plain Exception
                            try:
                                obj = getattr(obj, level)
                            except Exception:
                                break
                            if inspect.ismethod(obj):
                                obj = last_obj
                                class_attr, method = True, [level]
                                break
                    del remainder
                    is_class = inspect.isclass(obj)
                    if is_class or class_attr:
                        # Traverse all bases
                        classes = [obj if is_class else obj.__class__]
                        offset = 0
                        while offset < len(classes):
                            for base in classes[offset].__bases__:
                                # "object" as a base class is not very useful
                                if base not in classes and base is not object:
                                    classes.append(base)
                            offset += 1
                    else:
                        classes = [obj.__class__]
                    for cc in classes:
                        module = inspect.getmodule(cc)
                        if module is not None:
                            module = module.__name__.split(".")
                            class_name = cc.__qualname__
                            # a.b.C.meth could be documented as a.C.meth,
                            # so go down the list
                            for depth in range(len(module), 0, -1):
                                full_name = ".".join(
                                    module[:depth] + [class_name] + method
                                )
                                options.append(
                                    (name, full_name, class_attr, is_class, False)
                                )
            # second pass: by import (can't resolve as well without doing
            # some actions like actually importing the modules, so use it
            # as a last resort)
            for split_level in range(len(local_name_split)):
                local_name = ".".join(local_name_split[: split_level + 1])
                remainder = name[len(local_name) :]
                if local_name in self.imported_names:
                    full_name = self.imported_names[local_name] + remainder
                    is_class = class_attr = False  # can't tell without import
                    options.append((name, full_name, class_attr, is_class, False))
        return options


def _get_short_module_name(module_name, obj_name):
    """Get the shortest possible module name."""
    if "." in obj_name:
        obj_name, attr = obj_name.split(".")
    else:
        attr = None

    try:
        # look only in sys.modules to avoid importing the module, which may
        # otherwise have side effects
        real_obj = getattr(sys.modules[module_name], obj_name)
        if attr is not None:
            getattr(real_obj, attr)
    except (AttributeError, KeyError):
        # AttributeError: wrong class
        # KeyError: wrong object or module not previously imported
        return None

    parts = module_name.split(".")
    short_name = module_name
    for i in range(len(parts) - 1, 0, -1):
        short_name = ".".join(parts[:i])
        try:
            assert real_obj is getattr(sys.modules[short_name], obj_name)
        except (AssertionError, AttributeError, KeyError):
            # AssertionError: shortened object is not what we expect
            # KeyError: short module name not previously imported
            # AttributeError: wrong class or object
            # get the last working module name
            short_name = ".".join(parts[: (i + 1)])
            break
    return short_name


def _make_ref_regex(default_role=""):
    """Make regex to find reference to python objects."""
    # keep roles variable in sync values shown in configuration.rst
    # "Add mini-galleries for API documentation"
    roles = "func|meth|attr|obj|class"
    def_role_regex = (
        "|[^:]?" if re.fullmatch(f"(?:py:)?({roles})", default_role) else ""
    )
    # reference can have a separate title `title <reference>`,
    # don't match ``literal`` or `!disabled.references`
    return (
        rf"(?::(?:{roles}):{def_role_regex})"  # role prefix
        r"(?!``)`~?[^!]*?<?([^!~\s<>`]+)>?(?!``)`"
    )  # reference


def identify_names(script_blocks, ref_regex, global_variables=None, node=""):
    """Build a codeobj summary by identifying and resolving used names.

    Parameters
    ----------
    script_blocks : list
        (label, content, line_number)
        List where each element is a tuple with the label ('text' or 'code'),
        the corresponding content string of block and the leading line number.
    ref_regex : str
        Regex to find references to python objects.
    example_globals: Optional[Dict[str, Any]]
        Global variables for examples. Default=None
    node : ast.Module or str
        The parsed node. Default="".

    Returns
    -------
    example_code_obj : Dict[str, Any]
        Dict with information about all code object references found in an
        example. Dict contains the following keys:

            - example_code_obj['name'] : function or class name (str)
            - example_code_obj['module'] : module name (str)
            - example_code_obj['module_short'] : shortened module name (str)
            - example_code_obj['is_class'] : whether object is class (bool)
            - example_code_obj['is_explicit'] : whether object is an explicit
                backreference (referred to by sphinx markup) (bool)
    """
    if node == "":  # mostly convenience for testing functions
        c = "\n".join(block.content for block in script_blocks if block.type == "code")
        node = ast.parse(c)
    # Get matches from the code (AST, implicit matches)
    finder = NameFinder(global_variables)
    if node is not None:
        finder.visit(node)
    names = list(finder.get_mapping())
    # Get matches from docstring inspection (explicit matches)
    text = "\n".join(block.content for block in script_blocks if block.type == "text")
    names.extend((x, x, False, False, True) for x in re.findall(ref_regex, text))
    example_code_obj = dict()  # native dict preserves order nowadays
    # Make a list of all guesses, in `_embed_code_links` we will break
    # when we find a match
    for name, full_name, class_like, is_class, is_explicit in names:
        if name not in example_code_obj:
            example_code_obj[name] = list()
        # name is as written in file (e.g. np.asarray)
        # full_name includes resolved import path (e.g. numpy.asarray)
        splits = full_name.rsplit(".", 1 + class_like)
        if len(splits) == 1:
            splits = ("builtins", splits[0])
        elif len(splits) == 3:  # class-like
            assert class_like
            splits = (splits[0], ".".join(splits[1:]))
        else:
            assert not class_like

        module, attribute = splits

        # get shortened module name
        module_short = _get_short_module_name(module, attribute)
        cobj = dict(
            name=attribute,
            module=module,
            module_short=module_short or module,
            is_class=is_class,
            is_explicit=is_explicit,
        )
        example_code_obj[name].append(cobj)
    return example_code_obj


THUMBNAIL_TEMPLATE = """
.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="{intro}">

.. only:: html

  .. image:: /{thumbnail}
    :alt:

  :ref:`sphx_glr_{ref_name}`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">{title}</div>
    </div>

"""

BACKREF_THUMBNAIL_TEMPLATE = (
    THUMBNAIL_TEMPLATE
    + """
.. only:: not html

 * :ref:`sphx_glr_{ref_name}`
"""
)


def _thumbnail_div(
    target_dir, src_dir, fname, intro, title, is_backref=False, check=True
):
    """Generate reST to place a thumbnail in a gallery.

    Parameters
    ----------
    target_dir : str
        Absolute path to output directory, where thumbnails are saved.
    src_dir : str
        Absolute path to build source directory.
    fname : str
        Filename of example file.
    intro : str
        Introductory docstring of example, to show in tooltip
    title: str
        Title of example.
    is_backref : bool
        Whether the thumbnail is for a html backreference/recommendation file. If not
        a bullet list of the example is added for non-html documentation builds.
    check : bool
        Whether to check if the thumbnail image file exists.

    Returns
    -------
    thumbnail : str
        reST for a thumbnail.
    """
    fname = Path(fname)
    thumb, _ = _find_image_ext(
        os.path.join(target_dir, "images", "thumb", f"sphx_glr_{fname.stem}_thumb.png")
    )
    if check and not os.path.isfile(thumb):
        # This means we have done something wrong in creating our thumbnail!
        raise ExtensionError(
            f"Could not find internal Sphinx-Gallery thumbnail file:\n{thumb}"
        )
    thumb = os.path.relpath(thumb, src_dir)
    full_dir = os.path.relpath(target_dir, src_dir)

    # Inside rst files forward slash defines paths
    thumb = thumb.replace(os.sep, "/")

    ref_name = os.path.join(full_dir, fname).replace(os.sep, "_")

    template = BACKREF_THUMBNAIL_TEMPLATE if is_backref else THUMBNAIL_TEMPLATE
    return template.format(
        intro=escape(intro), thumbnail=thumb, title=title, ref_name=ref_name
    )


def _write_backreferences(
    backrefs, seen_backrefs, gallery_conf, src_dir, target_dir, fname, intro, title
):
    """Write and return backreferences for one example.

    Backreferences '.examples' file written includes reST of the list of examples
    as thumbnails.

    Parameters
    ----------
    backrefs : set[str]
        Back references to write.
    seen_backrefs: set
        Back references already encountered when parsing this example.
    gallery_conf : Dict[str, Any]
        Gallery configurations.
    src_dir : str
        Stuff.
    target_dir : str
        Absolute path to directory where examples are saved.
    fname : str
        Filename of current example python file.
    intro : str
        Introductory docstring of example.
    title: str
        Title of example.

    Returns
    -------
    backrefs_example : dict[str, tuple]
        Dictionary where value is the backreference object and value
        is a tuple containing: example filename, full path to example source directory,
        full path to example target directory, intro, title.
    """
    if gallery_conf["backreferences_dir"] is None:
        return

    backrefs_example = defaultdict(list)
    for backref in backrefs:
        include_path = os.path.join(
            gallery_conf["src_dir"],
            gallery_conf["backreferences_dir"],
            f"{backref}.examples.new",
        )
        seen = backref in seen_backrefs
        mode = "a" if seen else "w"
        with open(include_path, mode, **_W_KW) as ex_file:
            if not seen:
                # Be aware that if the number of lines of this heading changes,
                # the minigallery directive should be modified accordingly
                heading = f"Examples using ``{backref}``"
                ex_file.write("\n\n" + heading + "\n")
                ex_file.write("^" * len(heading) + "\n")
                ex_file.write("\n\n.. start-sphx-glr-thumbnails\n\n")
                # Open a div which will contain all thumbnails
                # (it will be closed in _finalize_backreferences)
                ex_file.write(THUMBNAIL_PARENT_DIV)
            ex_file.write(
                _thumbnail_div(
                    target_dir,
                    gallery_conf["src_dir"],
                    fname,
                    intro,
                    title,
                    is_backref=True,
                )
            )
            seen_backrefs.add(backref)
            backrefs_example[backref].append((fname, src_dir, target_dir, intro, title))
    return dict(backrefs_example)


def _finalize_backreferences(seen_backrefs, gallery_conf):
    """Replace backref files only if necessary."""
    logger = sphinx.util.logging.getLogger("sphinx-gallery")
    if gallery_conf["backreferences_dir"] is None:
        return

    for backref in seen_backrefs:
        path = os.path.join(
            gallery_conf["src_dir"],
            gallery_conf["backreferences_dir"],
            f"{backref}.examples.new",
        )
        if os.path.isfile(path):
            # Close div containing all thumbnails
            # (it was open in _write_backreferences)
            with open(path, "a", **_W_KW) as ex_file:
                ex_file.write(THUMBNAIL_PARENT_DIV_CLOSE)
            _replace_md5(path, mode="t")
        else:
            level = gallery_conf["log_level"]["backreference_missing"]
            func = getattr(logger, level)
            func(f"Could not find backreferences file: {path}")
            func(
                "The backreferences are likely to be erroneous "
                "due to file system case insensitivity."
            )
