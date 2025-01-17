"""Do syntax checks, but no writing."""

from __future__ import annotations

from typing import TYPE_CHECKING

from sphinx.builders import Builder
from sphinx.locale import __

if TYPE_CHECKING:
    from docutils import nodes

    from sphinx.application import Sphinx
    from sphinx.util.typing import ExtensionMetadata


class DummyBuilder(Builder):
    name = 'dummy'
    epilog = __('The dummy builder generates no files.')

    allow_parallel = True

    def init(self) -> None:
        pass

    def get_outdated_docs(self) -> set[str]:
        return self.env.found_docs

    def get_target_uri(self, docname: str, typ: str | None = None) -> str:
        return ''

    def write_doc(self, docname: str, doctree: nodes.document) -> None:
        pass

    def finish(self) -> None:
        pass


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_builder(DummyBuilder)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
