from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import os
import sys
from hashlib import md5

from docutils import nodes
from docutils.parsers.rst import directives
import warnings

from matplotlib import rcParams
from matplotlib.mathtext import MathTextParser
rcParams['mathtext.fontset'] = 'cm'
mathtext_parser = MathTextParser("Bitmap")

# Define LaTeX math node:
class latex_math(nodes.General, nodes.Element):
    pass

def fontset_choice(arg):
    return directives.choice(arg, ['cm', 'stix', 'stixsans'])

options_spec = {'fontset': fontset_choice}

def math_role(role, rawtext, text, lineno, inliner,
              options={}, content=[]):
    i = rawtext.find('`')
    latex = rawtext[i+1:-1]
    node = latex_math(rawtext)
    node['latex'] = latex
    node['fontset'] = options.get('fontset', 'cm')
    return [node], []
math_role.options = options_spec

def math_directive(name, arguments, options, content, lineno,
                   content_offset, block_text, state, state_machine):
    latex = ''.join(content)
    node = latex_math(block_text)
    node['latex'] = latex
    node['fontset'] = options.get('fontset', 'cm')
    return [node]

# This uses mathtext to render the expression
def latex2png(latex, filename, fontset='cm'):
    latex = "$%s$" % latex
    orig_fontset = rcParams['mathtext.fontset']
    rcParams['mathtext.fontset'] = fontset
    if os.path.exists(filename):
        depth = mathtext_parser.get_depth(latex, dpi=100)
    else:
        try:
            depth = mathtext_parser.to_png(filename, latex, dpi=100)
        except:
            warnings.warn("Could not render math expression %s" % latex,
                          Warning)
            depth = 0
    rcParams['mathtext.fontset'] = orig_fontset
    sys.stdout.write("#")
    sys.stdout.flush()
    return depth

# LaTeX to HTML translation stuff:
def latex2html(node, source):
    inline = isinstance(node.parent, nodes.TextElement)
    latex = node['latex']
    name = 'math-%s' % md5(latex.encode()).hexdigest()[-10:]

    destdir = os.path.join(setup.app.builder.outdir, '_images', 'mathmpl')
    if not os.path.exists(destdir):
        os.makedirs(destdir)
    dest = os.path.join(destdir, '%s.png' % name)
    path = '/'.join((setup.app.builder.imgpath, 'mathmpl'))

    depth = latex2png(latex, dest, node['fontset'])

    if inline:
        cls = ''
    else:
        cls = 'class="center" '
    if inline and depth != 0:
        style = 'style="position: relative; bottom: -%dpx"' % (depth + 1)
    else:
        style = ''

    return '<img src="%s/%s.png" %s%s/>' % (path, name, cls, style)


def setup(app):
    setup.app = app

    # Add visit/depart methods to HTML-Translator:
    def visit_latex_math_html(self, node):
        source = self.document.attributes['source']
        self.body.append(latex2html(node, source))

    def depart_latex_math_html(self, node):
        pass

    # Add visit/depart methods to LaTeX-Translator:
    def visit_latex_math_latex(self, node):
        inline = isinstance(node.parent, nodes.TextElement)
        if inline:
            self.body.append('$%s$' % node['latex'])
        else:
            self.body.extend(['\\begin{equation}',
                              node['latex'],
                              '\\end{equation}'])

    def depart_latex_math_latex(self, node):
        pass

    app.add_node(latex_math,
                 html=(visit_latex_math_html, depart_latex_math_html),
                 latex=(visit_latex_math_latex, depart_latex_math_latex))
    app.add_role('math', math_role)
    app.add_directive('math', math_directive,
                      True, (0, 0, 0), **options_spec)

    metadata = {'parallel_read_safe': True, 'parallel_write_safe': True}
    return metadata
