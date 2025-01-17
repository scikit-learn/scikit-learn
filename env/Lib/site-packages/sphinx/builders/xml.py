"""Docutils-native XML and pseudo-XML builders."""

from __future__ import annotations

from os import path
from typing import TYPE_CHECKING

from docutils import nodes
from docutils.io import StringOutput
from docutils.writers.docutils_xml import XMLTranslator

from sphinx.builders import Builder
from sphinx.locale import __
from sphinx.util import logging
from sphinx.util.osutil import (
    _last_modified_time,
    ensuredir,
    os_path,
)
from sphinx.writers.xml import PseudoXMLWriter, XMLWriter

if TYPE_CHECKING:
    from collections.abc import Iterator, Set

    from sphinx.application import Sphinx
    from sphinx.util.typing import ExtensionMetadata

logger = logging.getLogger(__name__)


class XMLBuilder(Builder):
    """
    Builds Docutils-native XML.
    """

    name = 'xml'
    format = 'xml'
    epilog = __('The XML files are in %(outdir)s.')

    out_suffix = '.xml'
    allow_parallel = True

    _writer_class: type[XMLWriter] | type[PseudoXMLWriter] = XMLWriter
    writer: XMLWriter | PseudoXMLWriter
    default_translator_class = XMLTranslator

    def init(self) -> None:
        pass

    def get_outdated_docs(self) -> Iterator[str]:
        for docname in self.env.found_docs:
            if docname not in self.env.all_docs:
                yield docname
                continue
            targetname = path.join(self.outdir, docname + self.out_suffix)
            try:
                targetmtime = _last_modified_time(targetname)
            except Exception:
                targetmtime = 0
            try:
                srcmtime = _last_modified_time(self.env.doc2path(docname))
                if srcmtime > targetmtime:
                    yield docname
            except OSError:
                # source doesn't exist anymore
                pass

    def get_target_uri(self, docname: str, typ: str | None = None) -> str:
        return docname

    def prepare_writing(self, docnames: Set[str]) -> None:
        self.writer = self._writer_class(self)

    def write_doc(self, docname: str, doctree: nodes.document) -> None:
        # work around multiple string % tuple issues in docutils;
        # replace tuples in attribute values with lists
        doctree = doctree.deepcopy()
        for domain in self.env.domains.sorted():
            doctree[f'xmlns:{domain.name}'] = 'https://www.sphinx-doc.org/'
        for node in doctree.findall(nodes.Element):
            for att, value in node.attributes.items():
                if isinstance(value, tuple):
                    node.attributes[att] = list(value)
                value = node.attributes[att]
                if isinstance(value, list):
                    for i, val in enumerate(value):
                        if isinstance(val, tuple):
                            value[i] = list(val)
        destination = StringOutput(encoding='utf-8')
        self.writer.write(doctree, destination)
        outfilename = path.join(self.outdir, os_path(docname) + self.out_suffix)
        ensuredir(path.dirname(outfilename))
        try:
            with open(outfilename, 'w', encoding='utf-8') as f:
                f.write(self.writer.output)
        except OSError as err:
            logger.warning(__('error writing file %s: %s'), outfilename, err)

    def finish(self) -> None:
        pass


class PseudoXMLBuilder(XMLBuilder):
    """
    Builds pseudo-XML for display purposes.
    """

    name = 'pseudoxml'
    format = 'pseudoxml'
    epilog = __('The pseudo-XML files are in %(outdir)s.')

    out_suffix = '.pseudoxml'

    _writer_class = PseudoXMLWriter


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_builder(XMLBuilder)
    app.add_builder(PseudoXMLBuilder)

    app.add_config_value('xml_pretty', True, 'env')

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
