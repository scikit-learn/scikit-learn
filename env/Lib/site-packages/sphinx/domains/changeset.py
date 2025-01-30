"""The changeset domain."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, ClassVar, NamedTuple

from docutils import nodes

from sphinx import addnodes
from sphinx.domains import Domain
from sphinx.locale import _
from sphinx.util.docutils import SphinxDirective

if TYPE_CHECKING:
    from collections.abc import Set

    from docutils.nodes import Node

    from sphinx.application import Sphinx
    from sphinx.environment import BuildEnvironment
    from sphinx.util.typing import ExtensionMetadata, OptionSpec


versionlabels = {
    'versionadded':   _('Added in version %s'),
    'versionchanged': _('Changed in version %s'),
    'deprecated':     _('Deprecated since version %s'),
    'versionremoved': _('Removed in version %s'),
}

versionlabel_classes = {
    'versionadded':     'added',
    'versionchanged':   'changed',
    'deprecated':       'deprecated',
    'versionremoved':   'removed',
}


class ChangeSet(NamedTuple):
    type: str
    docname: str
    lineno: int
    module: str | None
    descname: str | None
    content: str


class VersionChange(SphinxDirective):
    """
    Directive to describe a change/addition/deprecation in a specific version.
    """

    has_content = True
    required_arguments = 1
    optional_arguments = 1
    final_argument_whitespace = True
    option_spec: ClassVar[OptionSpec] = {}

    def run(self) -> list[Node]:
        node = addnodes.versionmodified()
        node.document = self.state.document
        self.set_source_info(node)
        node['type'] = self.name
        node['version'] = self.arguments[0]
        text = versionlabels[self.name] % self.arguments[0]
        if len(self.arguments) == 2:
            inodes, messages = self.parse_inline(self.arguments[1], lineno=self.lineno + 1)
            para = nodes.paragraph(self.arguments[1], '', *inodes, translatable=False)
            self.set_source_info(para)
            node.append(para)
        else:
            messages = []
        if self.content:
            node += self.parse_content_to_nodes()
        classes = ['versionmodified', versionlabel_classes[self.name]]
        if len(node) > 0 and isinstance(node[0], nodes.paragraph):
            # the contents start with a paragraph
            if node[0].rawsource:
                # make the first paragraph translatable
                content = nodes.inline(node[0].rawsource, translatable=True)
                content.source = node[0].source
                content.line = node[0].line
                content += node[0].children
                node[0].replace_self(nodes.paragraph('', '', content, translatable=False))

            para = node[0]
            para.insert(0, nodes.inline('', '%s: ' % text, classes=classes))
        elif len(node) > 0:
            # the contents do not starts with a paragraph
            para = nodes.paragraph('', '',
                                   nodes.inline('', '%s: ' % text, classes=classes),
                                   translatable=False)
            node.insert(0, para)
        else:
            # the contents are empty
            para = nodes.paragraph('', '',
                                   nodes.inline('', '%s.' % text, classes=classes),
                                   translatable=False)
            node.append(para)

        domain = self.env.domains.changeset_domain
        domain.note_changeset(node)

        ret: list[Node] = [node]
        ret += messages
        return ret


class ChangeSetDomain(Domain):
    """Domain for changesets."""

    name = 'changeset'
    label = 'changeset'

    initial_data: dict[str, dict[str, list[ChangeSet]]] = {
        'changes': {},      # version -> list of ChangeSet
    }

    @property
    def changesets(self) -> dict[str, list[ChangeSet]]:
        return self.data.setdefault('changes', {})  # version -> list of ChangeSet

    def note_changeset(self, node: addnodes.versionmodified) -> None:
        version = node['version']
        module = self.env.ref_context.get('py:module')
        objname = self.env.temp_data.get('object')
        changeset = ChangeSet(node['type'], self.env.docname, node.line,  # type: ignore[arg-type]
                              module, objname, node.astext())
        self.changesets.setdefault(version, []).append(changeset)

    def clear_doc(self, docname: str) -> None:
        for changes in self.changesets.values():
            for changeset in changes.copy():
                if changeset.docname == docname:
                    changes.remove(changeset)

    def merge_domaindata(self, docnames: Set[str], otherdata: dict[str, Any]) -> None:
        # XXX duplicates?
        for version, otherchanges in otherdata['changes'].items():
            changes = self.changesets.setdefault(version, [])
            for changeset in otherchanges:
                if changeset.docname in docnames:
                    changes.append(changeset)

    def process_doc(
        self, env: BuildEnvironment, docname: str, document: nodes.document,
    ) -> None:
        pass  # nothing to do here. All changesets are registered on calling directive.

    def get_changesets_for(self, version: str) -> list[ChangeSet]:
        return self.changesets.get(version, [])


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_domain(ChangeSetDomain)
    app.add_directive('deprecated', VersionChange)
    app.add_directive('versionadded', VersionChange)
    app.add_directive('versionchanged', VersionChange)
    app.add_directive('versionremoved', VersionChange)

    return {
        'version': 'builtin',
        'env_version': 1,
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
