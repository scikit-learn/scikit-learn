"""Domain indices."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, NamedTuple

from sphinx.errors import SphinxError

if TYPE_CHECKING:
    from collections.abc import Iterable

    from sphinx.domains import Domain


class IndexEntry(NamedTuple):
    """
    An index entry.

    .. note::

       The *qualifier* and *description* are not rendered for some output formats,
       such as LaTeX.
    """

    #: The name of the index entry to be displayed.
    name: str

    #: The sub-entry related type. One of:
    #:
    #: ``0``
    #:   A normal entry.
    #: ``1``
    #:   An entry with sub-entries.
    #: ``2``
    #:   A sub-entry.
    subtype: int

    #: *docname* where the entry is located.
    docname: str

    #: Anchor for the entry within `docname`
    anchor: str

    #: Extra info for the entry.
    extra: str

    #: Qualifier for the description.
    qualifier: str

    #: Description for the entry.
    descr: str


class Index(ABC):
    """
    An Index is the description for a domain-specific index.  To add an index to
    a domain, subclass Index, overriding the three name attributes:

    * `name` is an identifier used for generating file names.
      It is also used for a hyperlink target for the index. Therefore, users can
      refer the index page using ``ref`` role and a string which is combined
      domain name and ``name`` attribute (ex. ``:ref:`py-modindex```).
    * `localname` is the section title for the index.
    * `shortname` is a short name for the index, for use in the relation bar in
      HTML output.  Can be empty to disable entries in the relation bar.

    and providing a :meth:`generate()` method.  Then, add the index class to
    your domain's `indices` list.  Extensions can add indices to existing
    domains using :meth:`~sphinx.application.Sphinx.add_index_to_domain()`.

    .. versionchanged:: 3.0

       Index pages can be referred by domain name and index name via
       :rst:role:`ref` role.
    """

    name: str
    localname: str
    shortname: str | None = None

    def __init__(self, domain: Domain) -> None:
        if not self.name or self.localname is None:
            msg = f'Index subclass {self.__class__.__name__} has no valid name or localname'
            raise SphinxError(msg)
        self.domain = domain

    @abstractmethod
    def generate(
        self,
        docnames: Iterable[str] | None = None,
    ) -> tuple[list[tuple[str, list[IndexEntry]]], bool]:
        """Get entries for the index.

        If ``docnames`` is given, restrict to entries referring to these
        docnames.

        The return value is a tuple of ``(content, collapse)``:

        ``collapse``
          A boolean that determines if sub-entries should start collapsed (for
          output formats that support collapsing sub-entries).

        ``content``:
          A sequence of ``(letter, entries)`` tuples, where ``letter`` is the
          "heading" for the given ``entries``, usually the starting letter, and
          ``entries`` is a sequence of single entries.
          Each entry is an :py:class:`IndexEntry`.
        """
        raise NotImplementedError
