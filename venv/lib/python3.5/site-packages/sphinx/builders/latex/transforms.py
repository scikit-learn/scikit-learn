# -*- coding: utf-8 -*-
"""
    sphinx.builders.latex.transforms
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Transforms for LaTeX builder.

    :copyright: Copyright 2007-2018 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from docutils import nodes

from sphinx.transforms import SphinxTransform

if False:
    # For type annotation
    from typing import Dict, List, Set, Union # NOQA

URI_SCHEMES = ('mailto:', 'http:', 'https:', 'ftp:')


class ShowUrlsTransform(SphinxTransform):
    """Expand references to inline text or footnotes.

    For more information, see :confval:`latex_show_urls`.
    """
    default_priority = 400

    # references are expanded to footnotes (or not)
    expanded = False

    def apply(self):
        # type: () -> None
        # replace id_prefix temporarily
        id_prefix = self.document.settings.id_prefix
        self.document.settings.id_prefix = 'show_urls'

        self.expand_show_urls()
        if self.expanded:
            self.renumber_footnotes()

        # restore id_prefix
        self.document.settings.id_prefix = id_prefix

    def expand_show_urls(self):
        # type: () -> None
        show_urls = self.document.settings.env.config.latex_show_urls
        if show_urls is False or show_urls == 'no':
            return

        for node in self.document.traverse(nodes.reference):
            uri = node.get('refuri', '')
            if uri.startswith(URI_SCHEMES):
                if uri.startswith('mailto:'):
                    uri = uri[7:]
                if node.astext() != uri:
                    index = node.parent.index(node)
                    if show_urls == 'footnote':
                        footnote_nodes = self.create_footnote(uri)
                        for i, fn in enumerate(footnote_nodes):
                            node.parent.insert(index + i + 1, fn)

                        self.expanded = True
                    else:  # all other true values (b/w compat)
                        textnode = nodes.Text(" (%s)" % uri)
                        node.parent.insert(index + 1, textnode)

    def create_footnote(self, uri):
        # type: (unicode) -> List[Union[nodes.footnote, nodes.footnote_ref]]
        label = nodes.label('', '#')
        para = nodes.paragraph()
        para.append(nodes.reference('', nodes.Text(uri), refuri=uri, nolinkurl=True))
        footnote = nodes.footnote(uri, label, para, auto=1)
        footnote['names'].append('#')
        self.document.note_autofootnote(footnote)

        label = nodes.Text('#')
        footnote_ref = nodes.footnote_reference('[#]_', label, auto=1,
                                                refid=footnote['ids'][0])
        self.document.note_autofootnote_ref(footnote_ref)
        footnote.add_backref(footnote_ref['ids'][0])

        return [footnote, footnote_ref]

    def renumber_footnotes(self):
        # type: () -> None
        collector = FootnoteCollector(self.document)
        self.document.walkabout(collector)

        num = 0
        for document, footnote in collector.auto_footnotes:
            # search unused footnote number
            while True:
                num += 1
                if str(num) not in collector.used_footnote_numbers:
                    break

            # assign new footnote number
            old_label = footnote[0].astext()
            footnote[0].replace_self(nodes.label('', str(num)))
            if old_label in footnote['names']:
                footnote['names'].remove(old_label)
            footnote['names'].append(str(num))

            # update footnote_references by new footnote number
            for ref in collector.footnote_refs.get(document, []):
                if footnote['ids'][0] == ref['refid']:
                    ref.remove(ref[0])
                    ref += nodes.Text(str(num))


class FootnoteCollector(nodes.NodeVisitor):
    """Collect footnotes and footnote references on the document"""

    def __init__(self, document):
        # type: (nodes.document) -> None
        self.auto_footnotes = []            # type: List[nodes.footnote]
        self.used_footnote_numbers = set()  # type: Set[unicode]
        self.footnote_refs = {}             # type: Dict[nodes.Node, List[nodes.footnote_reference]]  # NOQA
        self.current_document = []          # type: List[nodes.Node]
        nodes.NodeVisitor.__init__(self, document)

    def visit_document(self, node):
        # type: (nodes.Node) -> None
        self.current_document.append(node)

    def depart_document(self, node):
        # type: (nodes.Node) -> None
        self.current_document.pop()

    def visit_start_of_file(self, node):
        # type: (nodes.Node) -> None
        self.current_document.append(node)

    def depart_start_of_file(self, node):
        # type: (nodes.Node) -> None
        self.current_document.pop()

    def unknown_visit(self, node):
        # type: (nodes.Node) -> None
        pass

    def visit_footnote(self, node):
        # type: (nodes.footnote) -> None
        document = self.current_document[-1]
        if node.get('auto'):
            self.auto_footnotes.append((document, node))
        else:
            for name in node['names']:
                self.used_footnote_numbers.add(name)

    def visit_footnote_reference(self, node):
        # type: (nodes.footnote_reference) -> None
        document = self.current_document[-1]
        footnote_refs = self.footnote_refs.setdefault(document, [])
        footnote_refs.append(node)

    def unknown_departure(self, node):
        # type: (nodes.Node) -> None
        pass
