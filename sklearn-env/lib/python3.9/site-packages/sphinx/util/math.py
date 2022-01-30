"""
    sphinx.util.math
    ~~~~~~~~~~~~~~~~

    Utility functions for math.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from docutils import nodes

from sphinx.builders.html import HTMLTranslator


def get_node_equation_number(writer: HTMLTranslator, node: nodes.math_block) -> str:
    if writer.builder.config.math_numfig and writer.builder.config.numfig:
        figtype = 'displaymath'
        if writer.builder.name == 'singlehtml':
            key = "%s/%s" % (writer.docnames[-1], figtype)
        else:
            key = figtype

        id = node['ids'][0]
        number = writer.builder.fignumbers.get(key, {}).get(id, ())
        return '.'.join(map(str, number))
    else:
        return node['number']


def wrap_displaymath(text: str, label: str, numbering: bool) -> str:
    def is_equation(part: str) -> str:
        return part.strip()

    if label is None:
        labeldef = ''
    else:
        labeldef = r'\label{%s}' % label
        numbering = True

    parts = list(filter(is_equation, text.split('\n\n')))
    equations = []
    if len(parts) == 0:
        return ''
    elif len(parts) == 1:
        if numbering:
            begin = r'\begin{equation}' + labeldef
            end = r'\end{equation}'
        else:
            begin = r'\begin{equation*}' + labeldef
            end = r'\end{equation*}'
        equations.append('\\begin{split}%s\\end{split}\n' % parts[0])
    else:
        if numbering:
            begin = r'\begin{align}%s\!\begin{aligned}' % labeldef
            end = r'\end{aligned}\end{align}'
        else:
            begin = r'\begin{align*}%s\!\begin{aligned}' % labeldef
            end = r'\end{aligned}\end{align*}'
        for part in parts:
            equations.append('%s\\\\\n' % part.strip())

    return '%s\n%s%s' % (begin, ''.join(equations), end)
