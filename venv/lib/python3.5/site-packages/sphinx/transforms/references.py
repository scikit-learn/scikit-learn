# -*- coding: utf-8 -*-
"""
    sphinx.transforms.references
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Docutils transforms used by Sphinx.

    :copyright: Copyright 2007-2018 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from docutils import nodes
from docutils.transforms.references import Substitutions

from sphinx.transforms import SphinxTransform


class SubstitutionDefinitionsRemover(SphinxTransform):
    """Remove ``substitution_definition node from doctrees.

    .. note:: In Sphinx-1.7, this transform is only used in LaTeX builder.
    """

    # should be invoked after Substitutions process
    default_priority = Substitutions.default_priority + 1

    def apply(self):
        # type: () -> None
        for node in self.document.traverse(nodes.substitution_definition):
            node.parent.remove(node)
