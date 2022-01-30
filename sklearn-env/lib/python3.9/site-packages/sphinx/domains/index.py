"""
    sphinx.domains.index
    ~~~~~~~~~~~~~~~~~~~~

    The index domain.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Tuple

from docutils import nodes
from docutils.nodes import Node, system_message
from docutils.parsers.rst import directives

from sphinx import addnodes
from sphinx.domains import Domain
from sphinx.environment import BuildEnvironment
from sphinx.util import logging, split_index_msg
from sphinx.util.docutils import ReferenceRole, SphinxDirective
from sphinx.util.nodes import process_index_entry
from sphinx.util.typing import OptionSpec

if TYPE_CHECKING:
    from sphinx.application import Sphinx


logger = logging.getLogger(__name__)


class IndexDomain(Domain):
    """Mathematics domain."""
    name = 'index'
    label = 'index'

    @property
    def entries(self) -> Dict[str, List[Tuple[str, str, str, str, str]]]:
        return self.data.setdefault('entries', {})

    def clear_doc(self, docname: str) -> None:
        self.entries.pop(docname, None)

    def merge_domaindata(self, docnames: Iterable[str], otherdata: Dict) -> None:
        for docname in docnames:
            self.entries[docname] = otherdata['entries'][docname]

    def process_doc(self, env: BuildEnvironment, docname: str, document: Node) -> None:
        """Process a document after it is read by the environment."""
        entries = self.entries.setdefault(env.docname, [])
        for node in list(document.findall(addnodes.index)):
            try:
                for entry in node['entries']:
                    split_index_msg(entry[0], entry[1])
            except ValueError as exc:
                logger.warning(str(exc), location=node)
                node.parent.remove(node)
            else:
                for entry in node['entries']:
                    entries.append(entry)


class IndexDirective(SphinxDirective):
    """
    Directive to add entries to the index.
    """
    has_content = False
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    option_spec: OptionSpec = {
        'name': directives.unchanged,
    }

    def run(self) -> List[Node]:
        arguments = self.arguments[0].split('\n')

        if 'name' in self.options:
            targetname = self.options['name']
            targetnode = nodes.target('', '', names=[targetname])
        else:
            targetid = 'index-%s' % self.env.new_serialno('index')
            targetnode = nodes.target('', '', ids=[targetid])

        self.state.document.note_explicit_target(targetnode)
        indexnode = addnodes.index()
        indexnode['entries'] = []
        indexnode['inline'] = False
        self.set_source_info(indexnode)
        for entry in arguments:
            indexnode['entries'].extend(process_index_entry(entry, targetnode['ids'][0]))
        return [indexnode, targetnode]


class IndexRole(ReferenceRole):
    def run(self) -> Tuple[List[Node], List[system_message]]:
        target_id = 'index-%s' % self.env.new_serialno('index')
        if self.has_explicit_title:
            # if an explicit target is given, process it as a full entry
            title = self.title
            entries = process_index_entry(self.target, target_id)
        else:
            # otherwise we just create a single entry
            if self.target.startswith('!'):
                title = self.title[1:]
                entries = [('single', self.target[1:], target_id, 'main', None)]
            else:
                title = self.title
                entries = [('single', self.target, target_id, '', None)]

        index = addnodes.index(entries=entries)
        target = nodes.target('', '', ids=[target_id])
        text = nodes.Text(title, title)
        self.set_source_info(index)
        return [index, target, text], []


def setup(app: "Sphinx") -> Dict[str, Any]:
    app.add_domain(IndexDomain)
    app.add_directive('index', IndexDirective)
    app.add_role('index', IndexRole())

    return {
        'version': 'builtin',
        'env_version': 1,
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
