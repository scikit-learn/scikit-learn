"""
    sphinx.ext.ifconfig
    ~~~~~~~~~~~~~~~~~~~

    Provides the ``ifconfig`` directive that allows to write documentation
    that is included depending on configuration variables.

    Usage::

        .. ifconfig:: releaselevel in ('alpha', 'beta', 'rc')

           This stuff is only included in the built docs for unstable versions.

    The argument for ``ifconfig`` is a plain Python expression, evaluated in the
    namespace of the project configuration (that is, all variables from
    ``conf.py`` are available.)

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from typing import Any, Dict, List

from docutils import nodes
from docutils.nodes import Node

import sphinx
from sphinx.application import Sphinx
from sphinx.util.docutils import SphinxDirective
from sphinx.util.nodes import nested_parse_with_titles
from sphinx.util.typing import OptionSpec


class ifconfig(nodes.Element):
    pass


class IfConfig(SphinxDirective):

    has_content = True
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    option_spec: OptionSpec = {}

    def run(self) -> List[Node]:
        node = ifconfig()
        node.document = self.state.document
        self.set_source_info(node)
        node['expr'] = self.arguments[0]
        nested_parse_with_titles(self.state, self.content, node)
        return [node]


def process_ifconfig_nodes(app: Sphinx, doctree: nodes.document, docname: str) -> None:
    ns = {confval.name: confval.value for confval in app.config}
    ns.update(app.config.__dict__.copy())
    ns['builder'] = app.builder.name
    for node in doctree.findall(ifconfig):
        try:
            res = eval(node['expr'], ns)
        except Exception as err:
            # handle exceptions in a clean fashion
            from traceback import format_exception_only
            msg = ''.join(format_exception_only(err.__class__, err))
            newnode = doctree.reporter.error('Exception occurred in '
                                             'ifconfig expression: \n%s' %
                                             msg, base_node=node)
            node.replace_self(newnode)
        else:
            if not res:
                node.replace_self([])
            else:
                node.replace_self(node.children)


def setup(app: Sphinx) -> Dict[str, Any]:
    app.add_node(ifconfig)
    app.add_directive('ifconfig', IfConfig)
    app.connect('doctree-resolved', process_ifconfig_nodes)
    return {'version': sphinx.__display_version__, 'parallel_read_safe': True}
