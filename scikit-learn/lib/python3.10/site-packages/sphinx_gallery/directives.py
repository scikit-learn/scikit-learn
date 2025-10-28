"""Custom Sphinx directives."""

import os
import shutil
from collections import namedtuple
from pathlib import Path, PurePosixPath

from docutils import nodes, statemachine
from docutils.parsers.rst import Directive, directives
from docutils.parsers.rst.directives import images
from sphinx.errors import ExtensionError
from sphinx.util.logging import getLogger

from .backreferences import (
    THUMBNAIL_PARENT_DIV,
    THUMBNAIL_PARENT_DIV_CLOSE,
    _thumbnail_div,
)
from .gen_rst import extract_intro_and_title
from .py_source_parser import split_code_and_text_blocks
from .utils import _read_json

logger = getLogger("sphinx-gallery")


class MiniGallery(Directive):
    """Custom directive to insert a mini-gallery.

    The required argument is one or more of the following:

    * fully qualified names of objects
    * pathlike strings to example Python files
    * glob-style pathlike strings to example Python files

    The string list of arguments is separated by spaces.

    The mini-gallery will be the subset of gallery
    examples that make use of that object from that specific namespace

    Options:

    * `add-heading` adds a heading to the mini-gallery.  If an argument is
      provided, it uses that text for the heading.  Otherwise, it uses
      default text.
    * `heading-level` specifies the heading level of the heading as a single
      character.  If omitted, the default heading level is `'^'`.
    """

    required_arguments = 0
    has_content = True
    optional_arguments = 1
    final_argument_whitespace = True
    option_spec = {
        "add-heading": directives.unchanged,
        "heading-level": directives.single_char_or_unicode,
    }

    def _get_target_dir(self, config, src_dir, path, obj):
        """Get thumbnail target directory, errors when not in example dir."""
        examples_dirs = config.sphinx_gallery_conf["examples_dirs"]
        if not isinstance(examples_dirs, list):
            examples_dirs = [examples_dirs]
        gallery_dirs = config.sphinx_gallery_conf["gallery_dirs"]
        if not isinstance(gallery_dirs, list):
            gallery_dirs = [gallery_dirs]

        gal_matches = []
        for e, g in zip(examples_dirs, gallery_dirs):
            e = Path(src_dir, e).resolve(strict=True)
            g = Path(src_dir, g).resolve(strict=True)
            try:
                gal_matches.append((path.relative_to(e), g))
            except ValueError:
                continue

        # `path` inside one `examples_dirs`
        if (n_match := len(gal_matches)) == 1:
            ex_parents, target_dir = (
                gal_matches[0][0].parents,
                gal_matches[0][1],
            )
        # `path` inside several `examples_dirs`, take the example dir with
        # longest match (shortest parents after relative)
        elif n_match > 1:
            ex_parents, target_dir = min(
                [(match[0].parents, match[1]) for match in gal_matches],
                key=lambda x: len(x[0]),
            )
        # `path` is not inside a `examples_dirs`
        else:
            raise ExtensionError(
                f"minigallery directive error: path input '{obj}' found file:"
                f" '{path}' but this does not live inside any examples_dirs: "
                f"{examples_dirs}"
            )

        # Add subgallery path, if present
        subdir = ""
        if (ex_p := ex_parents[0]) != Path("."):
            subdir = ex_p
        target_dir = target_dir / subdir
        return target_dir

    def run(self):
        """Generate mini-gallery from backreference and example files."""
        from .gen_rst import _get_callables

        if not (self.arguments or self.content):
            raise ExtensionError("No arguments passed to 'minigallery'")

        # Respect the same disabling options as the `raw` directive
        if (
            not self.state.document.settings.raw_enabled
            or not self.state.document.settings.file_insertion_enabled
        ):
            raise self.warning(f'"{self.name}" directive disabled.')

        # Retrieve the backreferences directory
        config = self.state.document.settings.env.config
        backreferences_dir = config.sphinx_gallery_conf["backreferences_dir"]
        if backreferences_dir is None:
            logger.warning(
                "'backreferences_dir' config is None, minigallery "
                "directive will resolve all inputs as file paths or globs."
            )

        # Retrieve source directory
        src_dir = config.sphinx_gallery_conf["src_dir"]

        # Parse the argument into the individual args
        arg_list = []

        if self.arguments:
            arg_list.extend([c.strip() for c in self.arguments[0].split()])
        if self.content:
            arg_list.extend([c.strip() for c in self.content])

        lines = []

        # Add a heading if requested
        if "add-heading" in self.options:
            heading = self.options["add-heading"]
            if heading == "":
                if len(arg_list) == 1:
                    heading = f"Examples using ``{arg_list[0]}``"
                else:
                    heading = "Examples using one of multiple objects"
            lines.append(heading)
            heading_level = self.options.get("heading-level", "^")
            lines.append(heading_level * len(heading))

        ExampleInfo = namedtuple("ExampleInfo", ["target_dir", "intro", "title", "arg"])
        backreferences_all = None
        if backreferences_dir:
            backreferences_all = _read_json(
                Path(src_dir, backreferences_dir, "backreferences_all.json")
            )

        def has_backrefs(arg):
            if backreferences_all is None:
                return False
            examples = backreferences_all.get(arg, None)
            if examples:
                return examples
            return False

        # Full paths to example files are dict keys, preventing duplicates
        file_paths = {}
        for arg in arg_list:
            # Backreference arg input
            if examples := has_backrefs(arg):
                # Note `ExampleInfo.arg` not required for backreference paths
                for path in examples:
                    file_paths[Path(path[1], path[0]).resolve()] = ExampleInfo(
                        target_dir=path[2], intro=path[3], title=path[4], arg=None
                    )
            # Glob path arg input
            elif paths := Path(src_dir).glob(arg):
                # Glob paths require extra parsing to get the intro and title
                # so we don't want to override a duplicate backreference arg input
                for path in paths:
                    path_resolved = path.resolve()
                    if path_resolved in file_paths:
                        continue
                    else:
                        file_paths[path_resolved] = ExampleInfo(None, None, None, arg)

        if len(file_paths) == 0:
            return []

        lines.append(THUMBNAIL_PARENT_DIV)

        # sort on the file path
        if config.sphinx_gallery_conf["minigallery_sort_order"] is None:
            sortkey = None
        else:
            (sortkey,) = _get_callables(
                config.sphinx_gallery_conf, "minigallery_sort_order"
            )
        for path, path_info in sorted(
            file_paths.items(),
            # `x[0]` to sort on key only
            key=((lambda x: sortkey(str(x[0]))) if sortkey else None),
        ):
            if path_info.intro is not None:
                thumbnail = _thumbnail_div(
                    path_info.target_dir,
                    src_dir,
                    path.name,
                    path_info.intro,
                    path_info.title,
                )
            else:
                target_dir = self._get_target_dir(config, src_dir, path, path_info.arg)
                # Get thumbnail
                # TODO: ideally we would not need to parse file (again) here
                _, script_blocks = split_code_and_text_blocks(
                    str(path), return_node=False
                )
                intro, title = extract_intro_and_title(
                    str(path), script_blocks[0].content
                )

                thumbnail = _thumbnail_div(target_dir, src_dir, path.name, intro, title)
            lines.append(thumbnail)

        lines.append(THUMBNAIL_PARENT_DIV_CLOSE)
        text = "\n".join(lines)
        include_lines = statemachine.string2lines(text, convert_whitespace=True)
        self.state_machine.insert_input(include_lines, str(path))

        return []


"""
Image sg for responsive images
"""


class imgsgnode(nodes.General, nodes.Element):
    """Sphinx Gallery image node class."""

    pass


class ImageSg(images.Image):
    """Implements a directive to allow an optional hidpi image.

    Meant to be used with the `image_srcset` configuration option.

    e.g.::

        .. image-sg:: /plot_types/basic/images/sphx_glr_bar_001.png
            :alt: bar
            :srcset: /plot_types/basic/images/sphx_glr_bar_001.png,
                     /plot_types/basic/images/sphx_glr_bar_001_2_00x.png 2.00x
            :class: sphx-glr-single-img

    The resulting html is::

        <img src="sphx_glr_bar_001_hidpi.png"
            srcset="_images/sphx_glr_bar_001.png,
                    _images/sphx_glr_bar_001_2_00x.png 2x",
            alt="bar"
            class="sphx-glr-single-img" />
    """

    has_content = False
    required_arguments = 1
    optional_arguments = 3
    final_argument_whitespace = False
    option_spec = {
        "srcset": directives.unchanged,
        "class": directives.class_option,
        "alt": directives.unchanged,
    }

    def run(self):
        """Update node contents."""
        image_node = imgsgnode()

        imagenm = self.arguments[0]
        image_node["alt"] = self.options.get("alt", "")
        image_node["class"] = self.options.get("class", None)

        # we would like uri to be the highest dpi version so that
        # latex etc will use that.  But for now, lets just make
        # imagenm

        image_node["uri"] = imagenm
        image_node["srcset"] = self.options.get("srcset", None)

        return [image_node]


def _parse_srcset(st):
    """Parse st."""
    entries = st.split(",")
    srcset = {}
    for entry in entries:
        spl = entry.strip().split(" ")
        if len(spl) == 1:
            srcset[0] = spl[0]
        elif len(spl) == 2:
            mult = spl[1][:-1]
            srcset[float(mult)] = spl[0]
        else:
            raise ExtensionError('srcset argument "{entry}" is invalid.')
    return srcset


def visit_imgsg_html(self, node):
    """Handle HTML image tag depending on 'srcset' configuration.

    If 'srcset' is not `None`, copy images, generate image html tag with 'srcset'
    and add to HTML `body`. If 'srcset' is `None` run `visit_image` on `node`.
    """
    if node["srcset"] is None:
        self.visit_image(node)
        return

    imagedir, srcset = _copy_images(self, node)

    # /doc/examples/subd/plot_1.rst
    docsource = self.document["source"]
    # /doc/
    # make sure to add the trailing slash:
    srctop = os.path.join(self.builder.srcdir, "")
    # examples/subd/plot_1.rst
    relsource = os.path.relpath(docsource, srctop)
    # /doc/build/html
    desttop = os.path.join(self.builder.outdir, "")
    # /doc/build/html/examples/subd
    dest = os.path.join(desttop, relsource)

    # ../../_images/ for dirhtml and ../_images/ for html
    imagerel = os.path.relpath(imagedir, os.path.dirname(dest))
    if self.builder.name == "dirhtml":
        imagerel = os.path.join("..", imagerel, "")
    else:  # html
        imagerel = os.path.join(imagerel, "")

    if "\\" in imagerel:
        imagerel = imagerel.replace("\\", "/")
    # make srcset str.  Need to change all the prefixes!
    srcsetst = ""
    for mult in srcset:
        nm = os.path.basename(srcset[mult][1:])
        # ../../_images/plot_1_2_0x.png
        relpath = imagerel + nm
        srcsetst += f"{relpath}"
        if mult == 0:
            srcsetst += ", "
        else:
            srcsetst += f" {mult:1.2f}x, "
    # trim trailing comma and space...
    srcsetst = srcsetst[:-2]

    # make uri also be relative...
    nm = os.path.basename(node["uri"][1:])
    uri = imagerel + nm

    alt = node["alt"]
    if node["class"] is not None:
        classst = node["class"][0]
        classst = f'class = "{classst}"'
    else:
        classst = ""

    html_block = f'<img src="{uri}" srcset="{srcsetst}" alt="{alt}"' + f" {classst}/>"
    self.body.append(html_block)


def visit_imgsg_latex(self, node):
    """Copy images, set node[uri] to highest resolution image and call `visit_image`."""
    if node["srcset"] is not None:
        imagedir, srcset = _copy_images(self, node)
        maxmult = -1
        # choose the highest res version for latex:
        for key in srcset.keys():
            maxmult = max(maxmult, key)
        node["uri"] = str(PurePosixPath(srcset[maxmult]).name)

    self.visit_image(node)


def _copy_images(self, node):
    srcset = _parse_srcset(node["srcset"])

    # where the sources are.  i.e. myproj/source
    srctop = self.builder.srcdir

    # copy image from source to imagedir.  This is
    # *probably* supposed to be done by a builder but...
    # ie myproj/build/html/_images
    imagedir = os.path.join(self.builder.imagedir, "")
    imagedir = PurePosixPath(self.builder.outdir, imagedir)

    os.makedirs(imagedir, exist_ok=True)

    # copy all the sources to the imagedir:
    for mult in srcset:
        abspath = PurePosixPath(srctop, srcset[mult][1:])
        shutil.copyfile(abspath, imagedir / abspath.name)

    return imagedir, srcset


def depart_imgsg_html(self, node):
    """HTML depart node visitor function."""
    pass


def depart_imgsg_latex(self, node):
    """LaTeX depart node visitor function."""
    self.depart_image(node)


def imagesg_addnode(app):
    """Add `imgsgnode` to Sphinx app with visitor functions for HTML and LaTeX."""
    app.add_node(
        imgsgnode,
        html=(visit_imgsg_html, depart_imgsg_html),
        latex=(visit_imgsg_latex, depart_imgsg_latex),
    )
