"""
    sphinx.builders.dummy
    ~~~~~~~~~~~~~~~~~~~~~

    Do syntax checks, but no writing.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from typing import Any, Dict, Set

from docutils.nodes import Node

from sphinx.application import Sphinx
from sphinx.builders import Builder
from sphinx.locale import __


class DummyBuilder(Builder):
    name = 'dummy'
    epilog = __('The dummy builder generates no files.')

    allow_parallel = True

    def init(self) -> None:
        pass

    def get_outdated_docs(self) -> Set[str]:
        return self.env.found_docs

    def get_target_uri(self, docname: str, typ: str = None) -> str:
        return ''

    def prepare_writing(self, docnames: Set[str]) -> None:
        pass

    def write_doc(self, docname: str, doctree: Node) -> None:
        pass

    def finish(self) -> None:
        pass


def setup(app: Sphinx) -> Dict[str, Any]:
    app.add_builder(DummyBuilder)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
