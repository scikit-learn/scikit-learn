"""
    sphinx.ext.autosectionlabel
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Allow reference sections by :ref: role using its title.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from typing import Any, Dict, cast

from docutils import nodes
from docutils.nodes import Node

from sphinx.application import Sphinx
from sphinx.domains.std import StandardDomain
from sphinx.locale import __
from sphinx.util import logging
from sphinx.util.nodes import clean_astext

logger = logging.getLogger(__name__)


def get_node_depth(node: Node) -> int:
    i = 0
    cur_node = node
    while cur_node.parent != node.document:
        cur_node = cur_node.parent
        i += 1
    return i


def register_sections_as_label(app: Sphinx, document: Node) -> None:
    domain = cast(StandardDomain, app.env.get_domain('std'))
    for node in document.findall(nodes.section):
        if (app.config.autosectionlabel_maxdepth and
                get_node_depth(node) >= app.config.autosectionlabel_maxdepth):
            continue
        labelid = node['ids'][0]
        docname = app.env.docname
        title = cast(nodes.title, node[0])
        ref_name = getattr(title, 'rawsource', title.astext())
        if app.config.autosectionlabel_prefix_document:
            name = nodes.fully_normalize_name(docname + ':' + ref_name)
        else:
            name = nodes.fully_normalize_name(ref_name)
        sectname = clean_astext(title)

        if name in domain.labels:
            logger.warning(__('duplicate label %s, other instance in %s'),
                           name, app.env.doc2path(domain.labels[name][0]),
                           location=node, type='autosectionlabel', subtype=docname)

        domain.anonlabels[name] = docname, labelid
        domain.labels[name] = docname, labelid, sectname


def setup(app: Sphinx) -> Dict[str, Any]:
    app.add_config_value('autosectionlabel_prefix_document', False, 'env')
    app.add_config_value('autosectionlabel_maxdepth', None, 'env')
    app.connect('doctree-read', register_sections_as_label)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
