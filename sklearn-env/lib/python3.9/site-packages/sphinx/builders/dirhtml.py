"""
    sphinx.builders.dirhtml
    ~~~~~~~~~~~~~~~~~~~~~~~

    Directory HTML builders.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from os import path
from typing import Any, Dict

from sphinx.application import Sphinx
from sphinx.builders.html import StandaloneHTMLBuilder
from sphinx.util import logging
from sphinx.util.osutil import SEP, os_path

logger = logging.getLogger(__name__)


class DirectoryHTMLBuilder(StandaloneHTMLBuilder):
    """
    A StandaloneHTMLBuilder that creates all HTML pages as "index.html" in
    a directory given by their pagename, so that generated URLs don't have
    ``.html`` in them.
    """
    name = 'dirhtml'

    def get_target_uri(self, docname: str, typ: str = None) -> str:
        if docname == 'index':
            return ''
        if docname.endswith(SEP + 'index'):
            return docname[:-5]  # up to sep
        return docname + SEP

    def get_outfilename(self, pagename: str) -> str:
        if pagename == 'index' or pagename.endswith(SEP + 'index'):
            outfilename = path.join(self.outdir, os_path(pagename) +
                                    self.out_suffix)
        else:
            outfilename = path.join(self.outdir, os_path(pagename),
                                    'index' + self.out_suffix)

        return outfilename


def setup(app: Sphinx) -> Dict[str, Any]:
    app.setup_extension('sphinx.builders.html')

    app.add_builder(DirectoryHTMLBuilder)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
