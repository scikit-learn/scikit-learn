"""
Custom Sphinx directives
========================
"""
import os
from pathlib import PurePosixPath
import shutil

from docutils import nodes
from docutils import statemachine
from docutils.parsers.rst import Directive, directives
from docutils.parsers.rst.directives import images

from sphinx.errors import ExtensionError


class MiniGallery(Directive):
    """
    Custom directive to insert a mini-gallery

    The required argument is one or more fully qualified names of objects,
    separated by spaces.  The mini-gallery will be the subset of gallery
    examples that make use of that object (from that specific namespace).

    Options:

    * `add-heading` adds a heading to the mini-gallery.  If an argument is
      provided, it uses that text for the heading.  Otherwise, it uses
      default text.
    * `heading-level` specifies the heading level of the heading as a single
      character.  If omitted, the default heading level is `'^'`.
    """
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    option_spec = {'add-heading': directives.unchanged,
                   'heading-level': directives.single_char_or_unicode}

    def run(self):
        # Respect the same disabling options as the `raw` directive
        if (not self.state.document.settings.raw_enabled
                or not self.state.document.settings.file_insertion_enabled):
            raise self.warning('"%s" directive disabled.' % self.name)

        # Retrieve the backreferences directory
        config = self.state.document.settings.env.config
        backreferences_dir = config.sphinx_gallery_conf['backreferences_dir']

        # Parse the argument into the individual objects
        obj_list = self.arguments[0].split()

        lines = []

        # Add a heading if requested
        if 'add-heading' in self.options:
            heading = self.options['add-heading']
            if heading == "":
                if len(obj_list) == 1:
                    heading = 'Examples using ``{}``'.format(obj_list[0])
                else:
                    heading = 'Examples using one of multiple objects'
            lines.append(heading)
            heading_level = self.options.get('heading-level', '^')
            lines.append(heading_level * len(heading))

        def has_backrefs(obj):
            src_dir = config.sphinx_gallery_conf['src_dir']
            path = os.path.join(src_dir, backreferences_dir,
                                '{}.examples'.format(obj))
            return os.path.isfile(path) and os.path.getsize(path) > 0

        if not any(has_backrefs(obj) for obj in obj_list):
            return []

        # Insert the backreferences file(s) using the `include` directive
        for obj in obj_list:
            path = os.path.join('/',  # Sphinx treats this as the source dir
                                backreferences_dir,
                                '{}.examples'.format(obj))

            # Always remove the heading (first 5 lines) from the file
            lines.append('.. include:: {}\n    :start-line: 5'.format(path))

        # Insert the end for the gallery using the `raw` directive
        lines.append('.. raw:: html\n\n    <div class="sphx-glr-clear"></div>')

        # Parse the assembly of `include` and `raw` directives
        text = '\n'.join(lines)
        include_lines = statemachine.string2lines(text,
                                                  convert_whitespace=True)
        self.state_machine.insert_input(include_lines, path)

        return []


"""
Image sg for responsive images
"""


class imgsgnode(nodes.General, nodes.Element):
    pass


def directive_boolean(value):
    if not value.strip():
        raise ValueError("No argument provided but required")
    if value.lower().strip() in ["yes", "1", 1, "true", "ok"]:
        return True
    elif value.lower().strip() in ['no', '0', 0, 'false', 'none']:
        return False
    else:
        raise ValueError(u"Please use one of: yes, true, no, false. "
                         u"Do not use `{}` as boolean.".format(value))


class ImageSg(images.Image):
    """
    Implements a directive to allow an optional hidpi image.  Meant to be
    used with the `image_srcset` configuration option.

    e.g.::

        .. image-sg:: /plot_types/basic/images/sphx_glr_bar_001.png
            :alt: bar
            :srcset: /plot_types/basic/images/sphx_glr_bar_001.png,
                     /plot_types/basic/images/sphx_glr_bar_001_2_0x.png 2.0x
            :class: sphx-glr-single-img

    The resulting html is::

        <img src="sphx_glr_bar_001_hidpi.png"
            srcset="_images/sphx_glr_bar_001.png,
                    _images/sphx_glr_bar_001_2_0x.png 2x",
            alt="bar"
            class="sphx-glr-single-img" />

    """

    has_content = False
    required_arguments = 1
    optional_arguments = 3
    final_argument_whitespace = False
    option_spec = {
        'srcset': directives.unchanged,
        'class': directives.class_option,
        'alt': directives.unchanged,
    }

    def run(self):
        image_node = imgsgnode()

        imagenm = self.arguments[0]
        image_node['alt'] = self.options.get('alt', '')
        image_node['class'] = self.options.get('class', None)

        # we would like uri to be the highest dpi version so that
        # latex etc will use that.  But for now, lets just make
        # imagenm

        image_node['uri'] = imagenm
        image_node['srcset'] = self.options.get('srcset', None)

        return [image_node]


def _parse_srcset(st):
    """ parse st"""
    entries = st.split(',')
    srcset = {}
    for entry in entries:
        spl = entry.strip().split(' ')
        if len(spl) == 1:
            srcset[0] = spl[0]
        elif len(spl) == 2:
            mult = spl[1][:-1]
            srcset[float(mult)] = spl[0]
        else:
            raise ExtensionError('srcset argument "{entry}" is invalid.')
    return srcset


def visit_imgsg_html(self, node):

    if node['srcset'] is None:
        self.visit_image(node)
        return

    imagedir, srcset = _copy_images(self, node)

    # /doc/examples/subd/plot_1.rst
    docsource = self.document['source']
    # /doc/
    # make sure to add the trailing slash:
    srctop = os.path.join(self.builder.srcdir, '')
    # examples/subd/plot_1.rst
    relsource = os.path.relpath(docsource, srctop)
    # /doc/build/html
    desttop = os.path.join(self.builder.outdir, '')
    # /doc/build/html/examples/subd
    dest = os.path.join(desttop, relsource)

    # ../../_images/
    imagerel = os.path.relpath(imagedir, os.path.dirname(dest))
    imagerel = os.path.join(imagerel, '')
    if '\\' in imagerel:
        imagerel = imagerel.replace('\\', '/')
    # make srcset str.  Need to change all the prefixes!
    srcsetst = ''
    for mult in srcset:
        nm = os.path.basename(srcset[mult][1:])
        # ../../_images/plot_1_2_0x.png
        relpath = imagerel+nm
        srcsetst += f'{relpath}'
        if mult == 0:
            srcsetst += ', '
        else:
            srcsetst += f' {mult:1.1f}x, '
    # trim trailing comma and space...
    srcsetst = srcsetst[:-2]

    # make uri also be relative...
    nm = os.path.basename(node['uri'][1:])
    uri = imagerel + nm

    alt = node['alt']
    if node['class'] is not None:
        classst = node['class'][0]
        classst = f'class = "{classst}"'
    else:
        classst = ''

    html_block = (f'<img src="{uri}" srcset="{srcsetst}" alt="{alt}"' +
                  f' {classst}/>')
    self.body.append(html_block)


def visit_imgsg_latex(self, node):

    if node['srcset'] is not None:
        imagedir, srcset = _copy_images(self, node)
        maxmult = -1
        # choose the highest res version for latex:
        for key in srcset.keys():
            maxmult = max(maxmult, key)
        node['uri'] = str(PurePosixPath(srcset[maxmult]).name)

    self.visit_image(node)


def _copy_images(self, node):
    srcset = _parse_srcset(node['srcset'])

    # where the sources are.  i.e. myproj/source
    srctop = self.builder.srcdir

    # copy image from source to imagedir.  This is
    # *probably* supposed to be done by a builder but...
    # ie myproj/build/html/_images
    imagedir = os.path.join(self.builder.imagedir, '')
    imagedir = PurePosixPath(self.builder.outdir, imagedir)

    os.makedirs(imagedir, exist_ok=True)

    # copy all the sources to the imagedir:
    for mult in srcset:
        abspath = PurePosixPath(srctop, srcset[mult][1:])
        shutil.copyfile(abspath, imagedir / abspath.name)

    return imagedir, srcset


def depart_imgsg_html(self, node):
    pass


def visit_sg_other(self, node):
    if node['uri'][0] == '/':
        node['uri'] = node['uri'][1:]
    self.visit_image(node)


def depart_imgsg_latex(self, node):
    self.depart_image(node)


def imagesg_addnode(app):
    app.add_node(imgsgnode,
                 html=(visit_imgsg_html, depart_imgsg_html),
                 latex=(visit_imgsg_latex, depart_imgsg_latex))
