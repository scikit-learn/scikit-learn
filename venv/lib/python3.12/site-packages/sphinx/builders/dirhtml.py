"""Directory HTML builders."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from sphinx.builders.html import StandaloneHTMLBuilder
from sphinx.util import logging
from sphinx.util.osutil import SEP

if TYPE_CHECKING:
    from sphinx.application import Sphinx
    from sphinx.util.typing import ExtensionMetadata

logger = logging.getLogger(__name__)


class DirectoryHTMLBuilder(StandaloneHTMLBuilder):
    """A StandaloneHTMLBuilder that creates all HTML pages as "index.html" in
    a directory given by their pagename, so that generated URLs don't have
    ``.html`` in them.
    """

    name = "dirhtml"

    def get_target_uri(self, docname: str, typ: str | None = None) -> str:
        if docname == "index":
            return ""
        if docname.endswith(SEP + "index"):
            return docname[:-5]  # up to sep
        return docname + SEP

    def get_output_path(self, page_name: str, /) -> Path:
        page_parts = page_name.split(SEP)
        if page_parts[-1] == "index":
            page_parts.pop()
        return Path(self.outdir, *page_parts, f"index{self.out_suffix}")


def setup(app: Sphinx) -> ExtensionMetadata:
    app.setup_extension("sphinx.builders.html")

    app.add_builder(DirectoryHTMLBuilder)

    return {
        "version": "builtin",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
