"""The dependencies collector components for sphinx.environment."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from sphinx.environment.collectors import EnvironmentCollector
from sphinx.util.osutil import _relative_path, fs_encoding

if TYPE_CHECKING:
    from collections.abc import Set

    from docutils import nodes

    from sphinx.application import Sphinx
    from sphinx.environment import BuildEnvironment
    from sphinx.util.typing import ExtensionMetadata


class DependenciesCollector(EnvironmentCollector):
    """dependencies collector for sphinx.environment."""

    def clear_doc(self, app: Sphinx, env: BuildEnvironment, docname: str) -> None:
        env.dependencies.pop(docname, None)

    def merge_other(
        self,
        app: Sphinx,
        env: BuildEnvironment,
        docnames: Set[str],
        other: BuildEnvironment,
    ) -> None:
        for docname in docnames:
            if docname in other.dependencies:
                env.dependencies[docname] = other.dependencies[docname]

    def process_doc(self, app: Sphinx, doctree: nodes.document) -> None:
        """Process docutils-generated dependency info."""
        cwd = Path.cwd()
        deps = doctree.settings.record_dependencies
        if not deps:
            return
        for dep in deps.list:
            # the dependency path is relative to the working dir, so get
            # one relative to the srcdir
            if isinstance(dep, bytes):
                dep = dep.decode(fs_encoding)
            relpath = _relative_path(cwd / dep, app.srcdir)
            app.env.note_dependency(relpath)


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_env_collector(DependenciesCollector)

    return {
        "version": "builtin",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
