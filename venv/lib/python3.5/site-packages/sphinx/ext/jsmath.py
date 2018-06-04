# -*- coding: utf-8 -*-
"""
    sphinx.ext.jsmath
    ~~~~~~~~~~~~~~~~~

    Set up everything for use of JSMath to display math in HTML
    via JavaScript.

    :copyright: Copyright 2007-2018 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from docutils import nodes

import sphinx
from sphinx.application import ExtensionError
from sphinx.ext.mathbase import get_node_equation_number
from sphinx.ext.mathbase import setup_math as mathbase_setup
from sphinx.locale import _


def html_visit_math(self, node):
    self.body.append(self.starttag(node, 'span', '', CLASS='math notranslate nohighlight'))
    self.body.append(self.encode(node['latex']) + '</span>')
    raise nodes.SkipNode


def html_visit_displaymath(self, node):
    if node['nowrap']:
        self.body.append(self.starttag(node, 'div', CLASS='math notranslate nohighlight'))
        self.body.append(self.encode(node['latex']))
        self.body.append('</div>')
        raise nodes.SkipNode
    for i, part in enumerate(node['latex'].split('\n\n')):
        part = self.encode(part)
        if i == 0:
            # necessary to e.g. set the id property correctly
            if node['number']:
                number = get_node_equation_number(self, node)
                self.body.append('<span class="eqno">(%s)' % number)
                self.add_permalink_ref(node, _('Permalink to this equation'))
                self.body.append('</span>')
            self.body.append(self.starttag(node, 'div', CLASS='math notranslate nohighlight'))
        else:
            # but only once!
            self.body.append('<div class="math">')
        if '&' in part or '\\\\' in part:
            self.body.append('\\begin{split}' + part + '\\end{split}')
        else:
            self.body.append(part)
        self.body.append('</div>\n')
    raise nodes.SkipNode


def builder_inited(app):
    if not app.config.jsmath_path:
        raise ExtensionError('jsmath_path config value must be set for the '
                             'jsmath extension to work')
    app.add_javascript(app.config.jsmath_path)


def setup(app):
    try:
        mathbase_setup(app, (html_visit_math, None), (html_visit_displaymath, None))
    except ExtensionError:
        raise ExtensionError('sphinx.ext.jsmath: other math package is already loaded')

    app.add_config_value('jsmath_path', '', False)
    app.connect('builder-inited', builder_inited)
    return {'version': sphinx.__display_version__, 'parallel_read_safe': True}
