"""
    sphinx.transforms.compact_bullet_list
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Docutils transforms used by Sphinx when reading documents.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from typing import Any, Dict, List, cast

from docutils import nodes
from docutils.nodes import Node

from sphinx import addnodes
from sphinx.application import Sphinx
from sphinx.transforms import SphinxTransform


class RefOnlyListChecker(nodes.GenericNodeVisitor):
    """Raise `nodes.NodeFound` if non-simple list item is encountered.

    Here 'simple' means a list item containing only a paragraph with a
    single reference in it.
    """

    def default_visit(self, node: Node) -> None:
        raise nodes.NodeFound

    def visit_bullet_list(self, node: nodes.bullet_list) -> None:
        pass

    def visit_list_item(self, node: nodes.list_item) -> None:
        children: List[Node] = []
        for child in node.children:
            if not isinstance(child, nodes.Invisible):
                children.append(child)
        if len(children) != 1:
            raise nodes.NodeFound
        if not isinstance(children[0], nodes.paragraph):
            raise nodes.NodeFound
        para = children[0]
        if len(para) != 1:
            raise nodes.NodeFound
        if not isinstance(para[0], addnodes.pending_xref):
            raise nodes.NodeFound
        raise nodes.SkipChildren

    def invisible_visit(self, node: Node) -> None:
        """Invisible nodes should be ignored."""
        pass


class RefOnlyBulletListTransform(SphinxTransform):
    """Change refonly bullet lists to use compact_paragraphs.

    Specifically implemented for 'Indices and Tables' section, which looks
    odd when html_compact_lists is false.
    """
    default_priority = 100

    def apply(self, **kwargs: Any) -> None:
        if self.config.html_compact_lists:
            return

        def check_refonly_list(node: Node) -> bool:
            """Check for list with only references in it."""
            visitor = RefOnlyListChecker(self.document)
            try:
                node.walk(visitor)
            except nodes.NodeFound:
                return False
            else:
                return True

        for node in self.document.findall(nodes.bullet_list):
            if check_refonly_list(node):
                for item in node.findall(nodes.list_item):
                    para = cast(nodes.paragraph, item[0])
                    ref = cast(nodes.reference, para[0])
                    compact_para = addnodes.compact_paragraph()
                    compact_para += ref
                    item.replace(para, compact_para)


def setup(app: Sphinx) -> Dict[str, Any]:
    app.add_transform(RefOnlyBulletListTransform)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
