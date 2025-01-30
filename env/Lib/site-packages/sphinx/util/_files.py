from __future__ import annotations

import hashlib
import os.path
from typing import Any


class FilenameUniqDict(dict[str, tuple[set[str], str]]):
    """
    A dictionary that automatically generates unique names for its keys,
    interpreted as filenames, and keeps track of a set of docnames they
    appear in.  Used for images and downloadable files in the environment.
    """

    def __init__(self) -> None:
        self._existing: set[str] = set()

    def add_file(self, docname: str, newfile: str) -> str:
        if newfile in self:
            self[newfile][0].add(docname)
            return self[newfile][1]
        uniquename = os.path.basename(newfile)
        base, ext = os.path.splitext(uniquename)
        i = 0
        while uniquename in self._existing:
            i += 1
            uniquename = f'{base}{i}{ext}'
        self[newfile] = ({docname}, uniquename)
        self._existing.add(uniquename)
        return uniquename

    def purge_doc(self, docname: str) -> None:
        for filename, (docs, unique) in list(self.items()):
            docs.discard(docname)
            if not docs:
                del self[filename]
                self._existing.discard(unique)

    def merge_other(
        self, docnames: set[str], other: dict[str, tuple[set[str], Any]]
    ) -> None:
        for filename, (docs, _unique) in other.items():
            for doc in docs & set(docnames):
                self.add_file(doc, filename)

    def __getstate__(self) -> set[str]:
        return self._existing

    def __setstate__(self, state: set[str]) -> None:
        self._existing = state


class DownloadFiles(dict[str, tuple[set[str], str]]):
    """A special dictionary for download files.

    .. important:: This class would be refactored in nearly future.
                   Hence don't hack this directly.
    """

    def add_file(self, docname: str, filename: str) -> str:
        if filename not in self:
            digest = hashlib.md5(filename.encode(), usedforsecurity=False).hexdigest()
            dest = f'{digest}/{os.path.basename(filename)}'
            self[filename] = (set(), dest)

        self[filename][0].add(docname)
        return self[filename][1]

    def purge_doc(self, docname: str) -> None:
        for filename, (docs, _dest) in list(self.items()):
            docs.discard(docname)
            if not docs:
                del self[filename]

    def merge_other(
        self, docnames: set[str], other: dict[str, tuple[set[str], Any]]
    ) -> None:
        for filename, (docs, _dest) in other.items():
            for docname in docs & set(docnames):
                self.add_file(docname, filename)
