"""
    sphinx.builders.html.transforms
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Transforms for HTML builder.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import re
from typing import Any, Dict, List

from docutils import nodes

from sphinx.application import Sphinx
from sphinx.transforms.post_transforms import SphinxPostTransform
from sphinx.util.nodes import NodeMatcher


class KeyboardTransform(SphinxPostTransform):
    """Transform :kbd: role to more detailed form.

    Before::

        <literal class="kbd">
            Control-x

    After::

        <literal class="kbd compound">
            <literal class="kbd">
                Control
            -
            <literal class="kbd">
                x
    """
    default_priority = 400
    formats = ('html',)
    pattern = re.compile(r'(?<=.)(-|\+|\^|\s+)(?=.)')
    multiwords_keys = (('caps', 'lock'),
                       ('page' 'down'),
                       ('page', 'up'),
                       ('scroll' 'lock'),
                       ('num', 'lock'),
                       ('sys' 'rq'),
                       ('back' 'space'))

    def run(self, **kwargs: Any) -> None:
        matcher = NodeMatcher(nodes.literal, classes=["kbd"])
        for node in self.document.findall(matcher):  # type: nodes.literal
            parts = self.pattern.split(node[-1].astext())
            if len(parts) == 1 or self.is_multiwords_key(parts):
                continue

            node['classes'].append('compound')
            node.pop()
            while parts:
                if self.is_multiwords_key(parts):
                    key = ''.join(parts[:3])
                    parts[:3] = []
                else:
                    key = parts.pop(0)
                node += nodes.literal('', key, classes=["kbd"])

                try:
                    # key separator (ex. -, +, ^)
                    sep = parts.pop(0)
                    node += nodes.Text(sep)
                except IndexError:
                    pass

    def is_multiwords_key(self, parts: List[str]) -> bool:
        if len(parts) >= 3 and parts[1].strip() == '':
            name = parts[0].lower(), parts[2].lower()
            if name in self.multiwords_keys:
                return True
            else:
                return False
        else:
            return False


def setup(app: Sphinx) -> Dict[str, Any]:
    app.add_post_transform(KeyboardTransform)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
