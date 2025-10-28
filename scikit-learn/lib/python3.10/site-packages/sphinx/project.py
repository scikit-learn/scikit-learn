"""Utility function and classes for Sphinx projects."""

from __future__ import annotations

import contextlib
import os
from pathlib import Path
from typing import TYPE_CHECKING

from sphinx.locale import __
from sphinx.util import logging
from sphinx.util._pathlib import _StrPath
from sphinx.util.matching import get_matching_files
from sphinx.util.osutil import path_stabilize

if TYPE_CHECKING:
    from collections.abc import Iterable

logger = logging.getLogger(__name__)
EXCLUDE_PATHS = ['**/_sources', '.#*', '**/.#*', '*.lproj/**']


class Project:
    """A project is the source code set of the Sphinx document(s)."""

    def __init__(
        self, srcdir: str | os.PathLike[str], source_suffix: Iterable[str]
    ) -> None:
        #: Source directory.
        self.srcdir = _StrPath(srcdir)

        #: source_suffix. Same as :confval:`source_suffix`.
        self.source_suffix = tuple(source_suffix)
        self._first_source_suffix = next(iter(self.source_suffix), '')

        #: The name of documents belonging to this project.
        self.docnames: set[str] = set()

        # Bijective mapping between docnames and (srcdir relative) paths.
        self._path_to_docname: dict[Path, str] = {}
        self._docname_to_path: dict[str, Path] = {}

    def restore(self, other: Project) -> None:
        """Take over a result of last build."""
        self.docnames = other.docnames
        self._path_to_docname = other._path_to_docname
        self._docname_to_path = other._docname_to_path

    def discover(
        self, exclude_paths: Iterable[str] = (), include_paths: Iterable[str] = ('**',)
    ) -> set[str]:
        """Find all document files in the source directory and put them in
        :attr:`docnames`.
        """
        self.docnames.clear()
        self._path_to_docname.clear()
        self._docname_to_path.clear()

        for filename in get_matching_files(
            self.srcdir,
            include_paths,
            [*exclude_paths, *EXCLUDE_PATHS],
        ):
            if docname := self.path2doc(filename):
                if docname in self.docnames:
                    files = [
                        str(f.relative_to(self.srcdir))
                        for f in self.srcdir.glob(f'{docname}.*')
                    ]
                    logger.warning(
                        __(
                            'multiple files found for the document "%s": %s\n'
                            'Use %r for the build.'
                        ),
                        docname,
                        ', '.join(files),
                        self.doc2path(docname, absolute=True),
                        once=True,
                    )
                elif os.access(self.srcdir / filename, os.R_OK):
                    self.docnames.add(docname)
                    path = Path(filename)
                    self._path_to_docname[path] = docname
                    self._docname_to_path[docname] = path
                else:
                    logger.warning(
                        __('Ignored unreadable document %r.'),
                        filename,
                        location=docname,
                    )

        return self.docnames

    def path2doc(self, filename: str | os.PathLike[str]) -> str | None:
        """Return the docname for the filename if the file is a document.

        *filename* should be absolute or relative to the source directory.
        """
        try:
            return self._path_to_docname[filename]  # type: ignore[index]
        except KeyError:
            path = Path(filename)
            if path.is_absolute():
                with contextlib.suppress(ValueError):
                    path = path.relative_to(self.srcdir)

            for suffix in self.source_suffix:
                if path.name.endswith(suffix):
                    return path_stabilize(path).removesuffix(suffix)

            # the file does not have a docname
            return None

    def doc2path(self, docname: str, absolute: bool) -> _StrPath:
        """Return the filename for the document name.

        If *absolute* is True, return as an absolute path.
        Else, return as a relative path to the source directory.
        """
        try:
            filename = self._docname_to_path[docname]
        except KeyError:
            # Backwards compatibility: the document does not exist
            filename = Path(docname + self._first_source_suffix)

        if absolute:
            return _StrPath(self.srcdir / filename)
        return _StrPath(filename)
