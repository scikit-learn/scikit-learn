"""
    sphinx.ext.graphviz
    ~~~~~~~~~~~~~~~~~~~

    Allow graphviz-formatted graphs to be included in Sphinx-generated
    documents inline.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import posixpath
import re
import subprocess
from os import path
from subprocess import PIPE, CalledProcessError
from typing import Any, Dict, List, Tuple

from docutils import nodes
from docutils.nodes import Node
from docutils.parsers.rst import Directive, directives

import sphinx
from sphinx.application import Sphinx
from sphinx.errors import SphinxError
from sphinx.locale import _, __
from sphinx.util import logging, sha1
from sphinx.util.docutils import SphinxDirective, SphinxTranslator
from sphinx.util.fileutil import copy_asset
from sphinx.util.i18n import search_image_for_language
from sphinx.util.nodes import set_source_info
from sphinx.util.osutil import ensuredir
from sphinx.util.typing import OptionSpec
from sphinx.writers.html import HTMLTranslator
from sphinx.writers.latex import LaTeXTranslator
from sphinx.writers.manpage import ManualPageTranslator
from sphinx.writers.texinfo import TexinfoTranslator
from sphinx.writers.text import TextTranslator

logger = logging.getLogger(__name__)


class GraphvizError(SphinxError):
    category = 'Graphviz error'


class ClickableMapDefinition:
    """A manipulator for clickable map file of graphviz."""
    maptag_re = re.compile('<map id="(.*?)"')
    href_re = re.compile('href=".*?"')

    def __init__(self, filename: str, content: str, dot: str = '') -> None:
        self.id: str = None
        self.filename = filename
        self.content = content.splitlines()
        self.clickable: List[str] = []

        self.parse(dot=dot)

    def parse(self, dot: str = None) -> None:
        matched = self.maptag_re.match(self.content[0])
        if not matched:
            raise GraphvizError('Invalid clickable map file found: %s' % self.filename)

        self.id = matched.group(1)
        if self.id == '%3':
            # graphviz generates wrong ID if graph name not specified
            # https://gitlab.com/graphviz/graphviz/issues/1327
            hashed = sha1(dot.encode()).hexdigest()
            self.id = 'grapviz%s' % hashed[-10:]
            self.content[0] = self.content[0].replace('%3', self.id)

        for line in self.content:
            if self.href_re.search(line):
                self.clickable.append(line)

    def generate_clickable_map(self) -> str:
        """Generate clickable map tags if clickable item exists.

        If not exists, this only returns empty string.
        """
        if self.clickable:
            return '\n'.join([self.content[0]] + self.clickable + [self.content[-1]])
        else:
            return ''


class graphviz(nodes.General, nodes.Inline, nodes.Element):
    pass


def figure_wrapper(directive: Directive, node: graphviz, caption: str) -> nodes.figure:
    figure_node = nodes.figure('', node)
    if 'align' in node:
        figure_node['align'] = node.attributes.pop('align')

    inodes, messages = directive.state.inline_text(caption, directive.lineno)
    caption_node = nodes.caption(caption, '', *inodes)
    caption_node.extend(messages)
    set_source_info(directive, caption_node)
    figure_node += caption_node
    return figure_node


def align_spec(argument: Any) -> str:
    return directives.choice(argument, ('left', 'center', 'right'))


class Graphviz(SphinxDirective):
    """
    Directive to insert arbitrary dot markup.
    """
    has_content = True
    required_arguments = 0
    optional_arguments = 1
    final_argument_whitespace = False
    option_spec: OptionSpec = {
        'alt': directives.unchanged,
        'align': align_spec,
        'caption': directives.unchanged,
        'layout': directives.unchanged,
        'graphviz_dot': directives.unchanged,  # an old alias of `layout` option
        'name': directives.unchanged,
        'class': directives.class_option,
    }

    def run(self) -> List[Node]:
        if self.arguments:
            document = self.state.document
            if self.content:
                return [document.reporter.warning(
                    __('Graphviz directive cannot have both content and '
                       'a filename argument'), line=self.lineno)]
            argument = search_image_for_language(self.arguments[0], self.env)
            rel_filename, filename = self.env.relfn2path(argument)
            self.env.note_dependency(rel_filename)
            try:
                with open(filename, encoding='utf-8') as fp:
                    dotcode = fp.read()
            except OSError:
                return [document.reporter.warning(
                    __('External Graphviz file %r not found or reading '
                       'it failed') % filename, line=self.lineno)]
        else:
            dotcode = '\n'.join(self.content)
            rel_filename = None
            if not dotcode.strip():
                return [self.state_machine.reporter.warning(
                    __('Ignoring "graphviz" directive without content.'),
                    line=self.lineno)]
        node = graphviz()
        node['code'] = dotcode
        node['options'] = {'docname': self.env.docname}

        if 'graphviz_dot' in self.options:
            node['options']['graphviz_dot'] = self.options['graphviz_dot']
        if 'layout' in self.options:
            node['options']['graphviz_dot'] = self.options['layout']
        if 'alt' in self.options:
            node['alt'] = self.options['alt']
        if 'align' in self.options:
            node['align'] = self.options['align']
        if 'class' in self.options:
            node['classes'] = self.options['class']
        if rel_filename:
            node['filename'] = rel_filename

        if 'caption' not in self.options:
            self.add_name(node)
            return [node]
        else:
            figure = figure_wrapper(self, node, self.options['caption'])
            self.add_name(figure)
            return [figure]


class GraphvizSimple(SphinxDirective):
    """
    Directive to insert arbitrary dot markup.
    """
    has_content = True
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec: OptionSpec = {
        'alt': directives.unchanged,
        'align': align_spec,
        'caption': directives.unchanged,
        'layout': directives.unchanged,
        'graphviz_dot': directives.unchanged,  # an old alias of `layout` option
        'name': directives.unchanged,
        'class': directives.class_option,
    }

    def run(self) -> List[Node]:
        node = graphviz()
        node['code'] = '%s %s {\n%s\n}\n' % \
                       (self.name, self.arguments[0], '\n'.join(self.content))
        node['options'] = {'docname': self.env.docname}
        if 'graphviz_dot' in self.options:
            node['options']['graphviz_dot'] = self.options['graphviz_dot']
        if 'layout' in self.options:
            node['options']['graphviz_dot'] = self.options['layout']
        if 'alt' in self.options:
            node['alt'] = self.options['alt']
        if 'align' in self.options:
            node['align'] = self.options['align']
        if 'class' in self.options:
            node['classes'] = self.options['class']

        if 'caption' not in self.options:
            self.add_name(node)
            return [node]
        else:
            figure = figure_wrapper(self, node, self.options['caption'])
            self.add_name(figure)
            return [figure]


def render_dot(self: SphinxTranslator, code: str, options: Dict, format: str,
               prefix: str = 'graphviz', filename: str = None) -> Tuple[str, str]:
    """Render graphviz code into a PNG or PDF output file."""
    graphviz_dot = options.get('graphviz_dot', self.builder.config.graphviz_dot)
    hashkey = (code + str(options) + str(graphviz_dot) +
               str(self.builder.config.graphviz_dot_args)).encode()

    fname = '%s-%s.%s' % (prefix, sha1(hashkey).hexdigest(), format)
    relfn = posixpath.join(self.builder.imgpath, fname)
    outfn = path.join(self.builder.outdir, self.builder.imagedir, fname)

    if path.isfile(outfn):
        return relfn, outfn

    if (hasattr(self.builder, '_graphviz_warned_dot') and
       self.builder._graphviz_warned_dot.get(graphviz_dot)):  # type: ignore  # NOQA
        return None, None

    ensuredir(path.dirname(outfn))

    dot_args = [graphviz_dot]
    dot_args.extend(self.builder.config.graphviz_dot_args)
    dot_args.extend(['-T' + format, '-o' + outfn])

    docname = options.get('docname', 'index')
    if filename:
        cwd = path.dirname(path.join(self.builder.srcdir, filename))
    else:
        cwd = path.dirname(path.join(self.builder.srcdir, docname))

    if format == 'png':
        dot_args.extend(['-Tcmapx', '-o%s.map' % outfn])

    try:
        ret = subprocess.run(dot_args, input=code.encode(), stdout=PIPE, stderr=PIPE,
                             cwd=cwd, check=True)
        if not path.isfile(outfn):
            raise GraphvizError(__('dot did not produce an output file:\n[stderr]\n%r\n'
                                   '[stdout]\n%r') % (ret.stderr, ret.stdout))
        return relfn, outfn
    except OSError:
        logger.warning(__('dot command %r cannot be run (needed for graphviz '
                          'output), check the graphviz_dot setting'), graphviz_dot)
        if not hasattr(self.builder, '_graphviz_warned_dot'):
            self.builder._graphviz_warned_dot = {}  # type: ignore
        self.builder._graphviz_warned_dot[graphviz_dot] = True  # type: ignore
        return None, None
    except CalledProcessError as exc:
        raise GraphvizError(__('dot exited with error:\n[stderr]\n%r\n'
                               '[stdout]\n%r') % (exc.stderr, exc.stdout)) from exc


def render_dot_html(self: HTMLTranslator, node: graphviz, code: str, options: Dict,
                    prefix: str = 'graphviz', imgcls: str = None, alt: str = None,
                    filename: str = None) -> Tuple[str, str]:
    format = self.builder.config.graphviz_output_format
    try:
        if format not in ('png', 'svg'):
            raise GraphvizError(__("graphviz_output_format must be one of 'png', "
                                   "'svg', but is %r") % format)
        fname, outfn = render_dot(self, code, options, format, prefix, filename)
    except GraphvizError as exc:
        logger.warning(__('dot code %r: %s'), code, exc)
        raise nodes.SkipNode from exc

    classes = [imgcls, 'graphviz'] + node.get('classes', [])
    imgcls = ' '.join(filter(None, classes))

    if fname is None:
        self.body.append(self.encode(code))
    else:
        if alt is None:
            alt = node.get('alt', self.encode(code).strip())
        if 'align' in node:
            self.body.append('<div align="%s" class="align-%s">' %
                             (node['align'], node['align']))
        if format == 'svg':
            self.body.append('<div class="graphviz">')
            self.body.append('<object data="%s" type="image/svg+xml" class="%s">\n' %
                             (fname, imgcls))
            self.body.append('<p class="warning">%s</p>' % alt)
            self.body.append('</object></div>\n')
        else:
            with open(outfn + '.map', encoding='utf-8') as mapfile:
                imgmap = ClickableMapDefinition(outfn + '.map', mapfile.read(), dot=code)
                if imgmap.clickable:
                    # has a map
                    self.body.append('<div class="graphviz">')
                    self.body.append('<img src="%s" alt="%s" usemap="#%s" class="%s" />' %
                                     (fname, alt, imgmap.id, imgcls))
                    self.body.append('</div>\n')
                    self.body.append(imgmap.generate_clickable_map())
                else:
                    # nothing in image map
                    self.body.append('<div class="graphviz">')
                    self.body.append('<img src="%s" alt="%s" class="%s" />' %
                                     (fname, alt, imgcls))
                    self.body.append('</div>\n')
        if 'align' in node:
            self.body.append('</div>\n')

    raise nodes.SkipNode


def html_visit_graphviz(self: HTMLTranslator, node: graphviz) -> None:
    render_dot_html(self, node, node['code'], node['options'], filename=node.get('filename'))


def render_dot_latex(self: LaTeXTranslator, node: graphviz, code: str,
                     options: Dict, prefix: str = 'graphviz', filename: str = None
                     ) -> None:
    try:
        fname, outfn = render_dot(self, code, options, 'pdf', prefix, filename)
    except GraphvizError as exc:
        logger.warning(__('dot code %r: %s'), code, exc)
        raise nodes.SkipNode from exc

    is_inline = self.is_inline(node)

    if not is_inline:
        pre = ''
        post = ''
        if 'align' in node:
            if node['align'] == 'left':
                pre = '{'
                post = r'\hspace*{\fill}}'
            elif node['align'] == 'right':
                pre = r'{\hspace*{\fill}'
                post = '}'
            elif node['align'] == 'center':
                pre = r'{\hfill'
                post = r'\hspace*{\fill}}'
        self.body.append('\n%s' % pre)

    self.body.append(r'\sphinxincludegraphics[]{%s}' % fname)

    if not is_inline:
        self.body.append('%s\n' % post)

    raise nodes.SkipNode


def latex_visit_graphviz(self: LaTeXTranslator, node: graphviz) -> None:
    render_dot_latex(self, node, node['code'], node['options'], filename=node.get('filename'))


def render_dot_texinfo(self: TexinfoTranslator, node: graphviz, code: str,
                       options: Dict, prefix: str = 'graphviz') -> None:
    try:
        fname, outfn = render_dot(self, code, options, 'png', prefix)
    except GraphvizError as exc:
        logger.warning(__('dot code %r: %s'), code, exc)
        raise nodes.SkipNode from exc
    if fname is not None:
        self.body.append('@image{%s,,,[graphviz],png}\n' % fname[:-4])
    raise nodes.SkipNode


def texinfo_visit_graphviz(self: TexinfoTranslator, node: graphviz) -> None:
    render_dot_texinfo(self, node, node['code'], node['options'])


def text_visit_graphviz(self: TextTranslator, node: graphviz) -> None:
    if 'alt' in node.attributes:
        self.add_text(_('[graph: %s]') % node['alt'])
    else:
        self.add_text(_('[graph]'))
    raise nodes.SkipNode


def man_visit_graphviz(self: ManualPageTranslator, node: graphviz) -> None:
    if 'alt' in node.attributes:
        self.body.append(_('[graph: %s]') % node['alt'])
    else:
        self.body.append(_('[graph]'))
    raise nodes.SkipNode


def on_build_finished(app: Sphinx, exc: Exception) -> None:
    if exc is None and app.builder.format == 'html':
        src = path.join(sphinx.package_dir, 'templates', 'graphviz', 'graphviz.css')
        dst = path.join(app.outdir, '_static')
        copy_asset(src, dst)


def setup(app: Sphinx) -> Dict[str, Any]:
    app.add_node(graphviz,
                 html=(html_visit_graphviz, None),
                 latex=(latex_visit_graphviz, None),
                 texinfo=(texinfo_visit_graphviz, None),
                 text=(text_visit_graphviz, None),
                 man=(man_visit_graphviz, None))
    app.add_directive('graphviz', Graphviz)
    app.add_directive('graph', GraphvizSimple)
    app.add_directive('digraph', GraphvizSimple)
    app.add_config_value('graphviz_dot', 'dot', 'html')
    app.add_config_value('graphviz_dot_args', [], 'html')
    app.add_config_value('graphviz_output_format', 'png', 'html')
    app.add_css_file('graphviz.css')
    app.connect('build-finished', on_build_finished)
    return {'version': sphinx.__display_version__, 'parallel_read_safe': True}
