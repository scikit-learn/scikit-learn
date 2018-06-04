#
# A pair of directives for inserting content that will only appear in
# either html or latex.
#

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from docutils.nodes import Body, Element


class only_base(Body, Element):
    def dont_traverse(self, *args, **kwargs):
        return []

class html_only(only_base):
    pass

class latex_only(only_base):
    pass

def run(content, node_class, state, content_offset):
    text = '\n'.join(content)
    node = node_class(text)
    state.nested_parse(content, content_offset, node)
    return [node]

def html_only_directive(name, arguments, options, content, lineno,
                        content_offset, block_text, state, state_machine):
    return run(content, html_only, state, content_offset)

def latex_only_directive(name, arguments, options, content, lineno,
                         content_offset, block_text, state, state_machine):
    return run(content, latex_only, state, content_offset)

def builder_inited(app):
    if app.builder.name == 'html':
        latex_only.traverse = only_base.dont_traverse
    else:
        html_only.traverse = only_base.dont_traverse


def setup(app):
    app.add_directive('htmlonly', html_only_directive, True, (0, 0, 0))
    app.add_directive('latexonly', latex_only_directive, True, (0, 0, 0))

    # This will *really* never see the light of day As it turns out,
    # this results in "broken" image nodes since they never get
    # processed, so best not to do this.
    # app.connect('builder-inited', builder_inited)

    # Add visit/depart methods to HTML-Translator:
    def visit_perform(self, node):
        pass

    def depart_perform(self, node):
        pass

    def visit_ignore(self, node):
        node.children = []

    def depart_ignore(self, node):
        node.children = []

    app.add_node(html_only,
                 html=(visit_perform, depart_perform),
                 latex=(visit_ignore, depart_ignore))
    app.add_node(latex_only,
                 latex=(visit_perform, depart_perform),
                 html=(visit_ignore, depart_ignore))

    metadata = {'parallel_read_safe': True, 'parallel_write_safe': True}
    return metadata
