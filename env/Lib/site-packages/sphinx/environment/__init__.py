"""Global creation environment."""

from __future__ import annotations

import functools
import os
import pickle
from collections import defaultdict
from copy import copy
from os import path
from typing import TYPE_CHECKING

from sphinx import addnodes
from sphinx.domains._domains_container import _DomainsContainer
from sphinx.environment.adapters import toctree as toctree_adapters
from sphinx.errors import (
    BuildEnvironmentError,
    DocumentError,
    ExtensionError,
    SphinxError,
)
from sphinx.locale import __
from sphinx.transforms import SphinxTransformer
from sphinx.util import logging
from sphinx.util._files import DownloadFiles, FilenameUniqDict
from sphinx.util._serialise import stable_str
from sphinx.util._timestamps import _format_rfc3339_microseconds
from sphinx.util.docutils import LoggingReporter
from sphinx.util.i18n import CatalogRepository, docname_to_domain
from sphinx.util.nodes import is_translatable
from sphinx.util.osutil import _last_modified_time, canon_path, os_path

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator
    from pathlib import Path
    from typing import Any, Literal

    from docutils import nodes
    from docutils.nodes import Node
    from docutils.parsers import Parser

    from sphinx.application import Sphinx
    from sphinx.builders import Builder
    from sphinx.config import Config
    from sphinx.domains import Domain
    from sphinx.events import EventManager
    from sphinx.project import Project
    from sphinx.util._pathlib import _StrPath

logger = logging.getLogger(__name__)

default_settings: dict[str, Any] = {
    'auto_id_prefix': 'id',
    'image_loading': 'link',
    'embed_stylesheet': False,
    'cloak_email_addresses': True,
    'pep_base_url': 'https://peps.python.org/',
    'pep_references': None,
    'rfc_base_url': 'https://datatracker.ietf.org/doc/html/',
    'rfc_references': None,
    'input_encoding': 'utf-8-sig',
    'doctitle_xform': False,
    'sectsubtitle_xform': False,
    'section_self_link': False,
    'halt_level': 5,
    'file_insertion_enabled': True,
    'smartquotes_locales': [],
}

# This is increased every time an environment attribute is added
# or changed to properly invalidate pickle files.
ENV_VERSION = 64

# config status
CONFIG_UNSET = -1
CONFIG_OK = 1
CONFIG_NEW = 2
CONFIG_CHANGED = 3
CONFIG_EXTENSIONS_CHANGED = 4

CONFIG_CHANGED_REASON = {
    CONFIG_NEW: __('new config'),
    CONFIG_CHANGED: __('config changed'),
    CONFIG_EXTENSIONS_CHANGED: __('extensions changed'),
}


versioning_conditions: dict[str, Literal[False] | Callable[[Node], bool]] = {
    'none': False,
    'text': is_translatable,
}


class BuildEnvironment:
    """
    The environment in which the ReST files are translated.
    Stores an inventory of cross-file targets and provides doctree
    transformations to resolve links to them.
    """

    # --------- ENVIRONMENT INITIALIZATION -------------------------------------

    def __init__(self, app: Sphinx) -> None:
        self.app: Sphinx = app
        self.doctreedir: Path = app.doctreedir
        self.srcdir: Path = app.srcdir
        self.config: Config = None  # type: ignore[assignment]
        self.config_status: int = CONFIG_UNSET
        self.config_status_extra: str = ''
        self.events: EventManager = app.events
        self.project: Project = app.project
        self.version: dict[str, int] = app.registry.get_envversion(app)

        # the method of doctree versioning; see set_versioning_method
        self.versioning_condition: Literal[False] | Callable[[Node], bool] | None = None
        self.versioning_compare: bool | None = None

        # the docutils settings for building
        self.settings: dict[str, Any] = default_settings.copy()
        self.settings['env'] = self

        # All "docnames" here are /-separated and relative and exclude
        # the source suffix.

        # docname -> time of reading (in integer microseconds)
        # contains all read docnames
        self.all_docs: dict[str, int] = {}
        # docname -> set of dependent file
        # names, relative to documentation root
        self.dependencies: dict[str, set[str]] = defaultdict(set)
        # docname -> set of included file
        # docnames included from other documents
        self.included: dict[str, set[str]] = defaultdict(set)
        # docnames to re-read unconditionally on next build
        self.reread_always: set[str] = set()

        self._pickled_doctree_cache: dict[str, bytes] = {}
        """In-memory cache for reading pickled doctrees from disk.
        docname -> pickled doctree

        This cache is used in the ``get_doctree`` method to avoid reading the
        doctree from disk multiple times.
        """

        self._write_doc_doctree_cache: dict[str, nodes.document] = {}
        """In-memory cache for unpickling doctrees from disk.
        docname -> doctree

        Items are added in ``Builder.write_doctree``, during the read phase,
        then used only in the ``get_and_resolve_doctree`` method.
        """

        # File metadata
        # docname -> dict of metadata items
        self.metadata: dict[str, dict[str, Any]] = defaultdict(dict)

        # TOC inventory
        # docname -> title node
        self.titles: dict[str, nodes.title] = {}
        # docname -> title node; only different if
        # set differently with title directive
        self.longtitles: dict[str, nodes.title] = {}
        # docname -> table of contents nodetree
        self.tocs: dict[str, nodes.bullet_list] = {}
        # docname -> number of real entries
        self.toc_num_entries: dict[str, int] = {}

        # used to determine when to show the TOC
        # in a sidebar (don't show if it's only one item)
        # docname -> dict of sectionid -> number
        self.toc_secnumbers: dict[str, dict[str, tuple[int, ...]]] = {}
        # docname -> dict of figtype -> dict of figureid -> number
        self.toc_fignumbers: dict[str, dict[str, dict[str, tuple[int, ...]]]] = {}

        # docname -> list of toctree includefiles
        self.toctree_includes: dict[str, list[str]] = {}
        # docname -> set of files (containing its TOCs) to rebuild too
        self.files_to_rebuild: dict[str, set[str]] = {}
        # docnames that have :glob: toctrees
        self.glob_toctrees: set[str] = set()
        # docnames that have :numbered: toctrees
        self.numbered_toctrees: set[str] = set()

        # domain-specific inventories, here to be pickled
        # domainname -> domain-specific dict
        self.domaindata: dict[str, dict[str, Any]] = {}

        # these map absolute path -> (docnames, unique filename)
        self.images: FilenameUniqDict = FilenameUniqDict()
        # filename -> (set of docnames, destination)
        self.dlfiles: DownloadFiles = DownloadFiles()

        # the original URI for images
        self.original_image_uri: dict[_StrPath, str] = {}

        # temporary data storage while reading a document
        self.temp_data: dict[str, Any] = {}
        # context for cross-references (e.g. current module or class)
        # this is similar to temp_data, but will for example be copied to
        # attributes of "any" cross references
        self.ref_context: dict[str, Any] = {}

        # search index data

        # docname -> title
        self._search_index_titles: dict[str, str | None] = {}
        # docname -> filename
        self._search_index_filenames: dict[str, str] = {}
        # stemmed words -> set(docname)
        self._search_index_mapping: dict[str, set[str]] = {}
        # stemmed words in titles -> set(docname)
        self._search_index_title_mapping: dict[str, set[str]] = {}
        # docname -> all titles in document
        self._search_index_all_titles: dict[str, list[tuple[str, str | None]]] = {}
        # docname -> list(index entry)
        self._search_index_index_entries: dict[str, list[tuple[str, str, str]]] = {}
        # objtype -> index
        self._search_index_objtypes: dict[tuple[str, str], int] = {}
        # objtype index -> (domain, type, objname (localized))
        self._search_index_objnames: dict[int, tuple[str, str, str]] = {}

        # all the registered domains, set by the application
        self.domains: _DomainsContainer = _DomainsContainer._from_environment(self)

        # set up environment
        self.setup(app)

    def __getstate__(self) -> dict[str, Any]:
        """Obtains serializable data for pickling."""
        __dict__ = self.__dict__.copy()
        # clear unpickable attributes
        __dict__.update(app=None, domains=None, events=None)
        # clear in-memory doctree caches, to reduce memory consumption and
        # ensure that, upon restoring the state, the most recent pickled files
        # on the disk are used instead of those from a possibly outdated state
        __dict__.update(_pickled_doctree_cache={}, _write_doc_doctree_cache={})
        return __dict__

    def __setstate__(self, state: dict[str, Any]) -> None:
        self.__dict__.update(state)

    def setup(self, app: Sphinx) -> None:
        """Set up BuildEnvironment object."""
        if self.version and self.version != app.registry.get_envversion(app):
            raise BuildEnvironmentError(__('build environment version not current'))
        if self.srcdir and self.srcdir != app.srcdir:
            raise BuildEnvironmentError(__('source directory has changed'))

        if self.project:
            app.project.restore(self.project)

        self.app = app
        self.doctreedir = app.doctreedir
        self.events = app.events
        self.srcdir = app.srcdir
        self.project = app.project
        self.version = app.registry.get_envversion(app)

        # initialise domains
        if self.domains is None:
            # if we are unpickling an environment, we need to recreate the domains
            self.domains = _DomainsContainer._from_environment(self)
        # setup domains (must do after all initialization)
        self.domains._setup()

        # Initialise config.
        # The old config is self.config, restored from the pickled environment.
        # The new config is app.config, always recreated from ``conf.py``
        self.config_status, self.config_status_extra = self._config_status(
            old_config=self.config, new_config=app.config, verbosity=app.verbosity
        )
        self.config = app.config

        # initialize settings
        self._update_settings(app.config)

    @staticmethod
    def _config_status(
        *, old_config: Config | None, new_config: Config, verbosity: int
    ) -> tuple[int, str]:
        """Report the differences between two Config objects.

        Returns a triple of:

        1. The new configuration
        2. A status code indicating how the configuration has changed.
        3. A status message indicating what has changed.
        """
        if old_config is None:
            return CONFIG_NEW, ''

        if old_config.extensions != new_config.extensions:
            old_extensions = set(old_config.extensions)
            new_extensions = set(new_config.extensions)
            extensions = old_extensions ^ new_extensions
            if len(extensions) == 1:
                extension = extensions.pop()
            else:
                extension = f'{len(extensions)}'
            return CONFIG_EXTENSIONS_CHANGED, f' ({extension!r})'

        # Log any changes in configuration keys
        if changed_keys := _differing_config_keys(old_config, new_config):
            changed_num = len(changed_keys)
            if changed_num == 1:
                logger.info(
                    __('The configuration has changed (1 option: %r)'),
                    next(iter(changed_keys)),
                )
            elif changed_num <= 5 or verbosity >= 1:
                logger.info(
                    __('The configuration has changed (%d options: %s)'),
                    changed_num,
                    ', '.join(map(repr, sorted(changed_keys))),
                )
            else:
                logger.info(
                    __('The configuration has changed (%d options: %s, ...)'),
                    changed_num,
                    ', '.join(map(repr, sorted(changed_keys)[:5])),
                )

        # check if a config value was changed that affects how doctrees are read
        for item in new_config.filter(frozenset({'env'})):
            if old_config[item.name] != item.value:
                return CONFIG_CHANGED, f' ({item.name!r})'

        return CONFIG_OK, ''

    def _update_settings(self, config: Config) -> None:
        """Update settings by new config."""
        self.settings['input_encoding'] = config.source_encoding
        self.settings['trim_footnote_reference_space'] = (
            config.trim_footnote_reference_space
        )
        self.settings['language_code'] = config.language

        # Allow to disable by 3rd party extension (workaround)
        self.settings.setdefault('smart_quotes', True)

    def set_versioning_method(
        self, method: str | Callable[[Node], bool], compare: bool
    ) -> None:
        """Set the doctree versioning method for this environment.

        Versioning methods are a builder property; only builders with the same
        versioning method can share the same doctree directory.  Therefore, we
        raise an exception if the user tries to use an environment with an
        incompatible versioning method.
        """
        condition: Literal[False] | Callable[[Node], bool]
        if callable(method):
            condition = method
        else:
            if method not in versioning_conditions:
                raise ValueError('invalid versioning method: %r' % method)
            condition = versioning_conditions[method]

        if self.versioning_condition not in {None, condition}:
            msg = __(
                'This environment is incompatible with the '
                'selected builder, please choose another '
                'doctree directory.'
            )
            raise SphinxError(msg)
        self.versioning_condition = condition
        self.versioning_compare = compare

    def clear_doc(self, docname: str) -> None:
        """Remove all traces of a source file in the inventory."""
        if docname in self.all_docs:
            self.all_docs.pop(docname, None)
            self.included.pop(docname, None)
            self.reread_always.discard(docname)

        self.domains._clear_doc(docname)

    def merge_info_from(
        self, docnames: Iterable[str], other: BuildEnvironment, app: Sphinx
    ) -> None:
        """Merge global information gathered about *docnames* while reading them
        from the *other* environment.

        This possibly comes from a parallel build process.
        """
        docnames = frozenset(docnames)
        for docname in docnames:
            self.all_docs[docname] = other.all_docs[docname]
            self.included[docname] = other.included[docname]
            if docname in other.reread_always:
                self.reread_always.add(docname)

        self.domains._merge_domain_data(docnames, other.domaindata)
        self.events.emit('env-merge-info', self, docnames, other)

    def path2doc(self, filename: str | os.PathLike[str]) -> str | None:
        """Return the docname for the filename if the file is document.

        *filename* should be absolute or relative to the source directory.
        """
        return self.project.path2doc(filename)

    def doc2path(self, docname: str, base: bool = True) -> _StrPath:
        """Return the filename for the document name.

        If *base* is True, return absolute path under self.srcdir.
        If *base* is False, return relative path to self.srcdir.
        """
        return self.project.doc2path(docname, absolute=base)

    def relfn2path(self, filename: str, docname: str | None = None) -> tuple[str, str]:
        """Return paths to a file referenced from a document, relative to
        documentation root and absolute.

        In the input "filename", absolute filenames are taken as relative to the
        source dir, while relative filenames are relative to the dir of the
        containing document.
        """
        filename = os_path(filename)
        if filename.startswith(('/', os.sep)):
            rel_fn = filename[1:]
        else:
            docdir = path.dirname(self.doc2path(docname or self.docname, base=False))
            rel_fn = path.join(docdir, filename)

        return (
            canon_path(path.normpath(rel_fn)),
            path.normpath(path.join(self.srcdir, rel_fn)),
        )

    @property
    def found_docs(self) -> set[str]:
        """Contains all existing docnames."""
        return self.project.docnames

    def find_files(self, config: Config, builder: Builder) -> None:
        """Find all source files in the source dir and put them in
        self.found_docs.
        """
        try:
            exclude_paths = (
                self.config.exclude_patterns
                + self.config.templates_path
                + builder.get_asset_paths()
            )
            self.project.discover(exclude_paths, self.config.include_patterns)

            # Current implementation is applying translated messages in the reading
            # phase.Therefore, in order to apply the updated message catalog, it is
            # necessary to re-process from the reading phase. Here, if dependency
            # is set for the doc source and the mo file, it is processed again from
            # the reading phase when mo is updated. In the future, we would like to
            # move i18n process into the writing phase, and remove these lines.
            if builder.use_message_catalog:
                # add catalog mo file dependency
                repo = CatalogRepository(
                    self.srcdir,
                    self.config.locale_dirs,
                    self.config.language,
                    self.config.source_encoding,
                )
                mo_paths = {c.domain: c.mo_path for c in repo.catalogs}
                for docname in self.found_docs:
                    domain = docname_to_domain(docname, self.config.gettext_compact)
                    if domain in mo_paths:
                        self.dependencies[docname].add(mo_paths[domain])
        except OSError as exc:
            raise DocumentError(
                __('Failed to scan documents in %s: %r') % (self.srcdir, exc)
            ) from exc

    def get_outdated_files(
        self, config_changed: bool
    ) -> tuple[set[str], set[str], set[str]]:
        """Return (added, changed, removed) sets."""
        # clear all files no longer present
        removed = set(self.all_docs) - self.found_docs

        added: set[str] = set()
        changed: set[str] = set()

        if config_changed:
            # config values affect e.g. substitutions
            added = self.found_docs
        else:
            for docname in self.found_docs:
                if docname not in self.all_docs:
                    logger.debug('[build target] added %r', docname)
                    added.add(docname)
                    continue
                # if the doctree file is not there, rebuild
                filename = path.join(self.doctreedir, docname + '.doctree')
                if not path.isfile(filename):
                    logger.debug('[build target] changed %r', docname)
                    changed.add(docname)
                    continue
                # check the "reread always" list
                if docname in self.reread_always:
                    logger.debug('[build target] changed %r', docname)
                    changed.add(docname)
                    continue
                # check the mtime of the document
                mtime = self.all_docs[docname]
                newmtime = _last_modified_time(self.doc2path(docname))
                if newmtime > mtime:
                    logger.debug(
                        '[build target] outdated %r: %s -> %s',
                        docname,
                        _format_rfc3339_microseconds(mtime),
                        _format_rfc3339_microseconds(newmtime),
                    )
                    changed.add(docname)
                    continue
                # finally, check the mtime of dependencies
                for dep in self.dependencies[docname]:
                    try:
                        # this will do the right thing when dep is absolute too
                        deppath = path.join(self.srcdir, dep)
                        if not path.isfile(deppath):
                            logger.debug(
                                '[build target] changed %r missing dependency %r',
                                docname,
                                deppath,
                            )
                            changed.add(docname)
                            break
                        depmtime = _last_modified_time(deppath)
                        if depmtime > mtime:
                            logger.debug(
                                '[build target] outdated %r from dependency %r: %s -> %s',
                                docname,
                                deppath,
                                _format_rfc3339_microseconds(mtime),
                                _format_rfc3339_microseconds(depmtime),
                            )
                            changed.add(docname)
                            break
                    except OSError:
                        # give it another chance
                        changed.add(docname)
                        break

        return added, changed, removed

    def check_dependents(self, app: Sphinx, already: set[str]) -> Iterator[str]:
        to_rewrite: list[str] = []
        for docnames in self.events.emit('env-get-updated', self):
            to_rewrite.extend(docnames)
        for docname in set(to_rewrite):
            if docname not in already:
                yield docname

    # --------- SINGLE FILE READING --------------------------------------------

    def prepare_settings(self, docname: str) -> None:
        """Prepare to set up environment for reading."""
        self.temp_data['docname'] = docname
        # defaults to the global default, but can be re-set in a document
        self.temp_data['default_role'] = self.config.default_role
        self.temp_data['default_domain'] = self.domains.get(self.config.primary_domain)

    # utilities to use while reading a document

    @property
    def docname(self) -> str:
        """Returns the docname of the document currently being parsed."""
        return self.temp_data['docname']

    @property
    def parser(self) -> Parser:
        """Returns the parser being used for to parse the current document."""
        return self.temp_data['_parser']

    def new_serialno(self, category: str = '') -> int:
        """Return a serial number, e.g. for index entry targets.

        The number is guaranteed to be unique in the current document.
        """
        key = category + 'serialno'
        cur = self.temp_data.get(key, 0)
        self.temp_data[key] = cur + 1
        return cur

    def note_dependency(self, filename: str) -> None:
        """Add *filename* as a dependency of the current document.

        This means that the document will be rebuilt if this file changes.

        *filename* should be absolute or relative to the source directory.
        """
        self.dependencies[self.docname].add(filename)

    def note_included(self, filename: str) -> None:
        """Add *filename* as a included from other document.

        This means the document is not orphaned.

        *filename* should be absolute or relative to the source directory.
        """
        doc = self.path2doc(filename)
        if doc:
            self.included[self.docname].add(doc)

    def note_reread(self) -> None:
        """Add the current document to the list of documents that will
        automatically be re-read at the next build.
        """
        self.reread_always.add(self.docname)

    def get_domain(self, domainname: str) -> Domain:
        """Return the domain instance with the specified name.

        Raises an ExtensionError if the domain is not registered.
        """
        try:
            return self.domains[domainname]
        except KeyError as exc:
            msg = __('Domain %r is not registered') % domainname
            raise ExtensionError(msg) from exc

    # --------- RESOLVING REFERENCES AND TOCTREES ------------------------------

    def get_doctree(self, docname: str) -> nodes.document:
        """Read the doctree for a file from the pickle and return it."""
        try:
            serialised = self._pickled_doctree_cache[docname]
        except KeyError:
            filename = path.join(self.doctreedir, docname + '.doctree')
            with open(filename, 'rb') as f:
                serialised = self._pickled_doctree_cache[docname] = f.read()

        doctree = pickle.loads(serialised)
        doctree.settings.env = self
        doctree.reporter = LoggingReporter(str(self.doc2path(docname)))
        return doctree

    @functools.cached_property
    def master_doctree(self) -> nodes.document:
        return self.get_doctree(self.config.root_doc)

    def get_and_resolve_doctree(
        self,
        docname: str,
        builder: Builder,
        doctree: nodes.document | None = None,
        prune_toctrees: bool = True,
        includehidden: bool = False,
    ) -> nodes.document:
        """Read the doctree from the pickle, resolve cross-references and
        toctrees and return it.
        """
        if doctree is None:
            try:
                doctree = self._write_doc_doctree_cache.pop(docname)
                doctree.settings.env = self
                doctree.reporter = LoggingReporter(str(self.doc2path(docname)))
            except KeyError:
                doctree = self.get_doctree(docname)

        # resolve all pending cross-references
        self.apply_post_transforms(doctree, docname)

        # now, resolve all toctree nodes
        for toctreenode in doctree.findall(addnodes.toctree):
            result = toctree_adapters._resolve_toctree(
                self,
                docname,
                builder,
                toctreenode,
                prune=prune_toctrees,
                includehidden=includehidden,
            )
            if result is None:
                toctreenode.parent.replace(toctreenode, [])
            else:
                toctreenode.replace_self(result)

        return doctree

    def resolve_toctree(
        self,
        docname: str,
        builder: Builder,
        toctree: addnodes.toctree,
        prune: bool = True,
        maxdepth: int = 0,
        titles_only: bool = False,
        collapse: bool = False,
        includehidden: bool = False,
    ) -> Node | None:
        """Resolve a *toctree* node into individual bullet lists with titles
        as items, returning None (if no containing titles are found) or
        a new node.

        If *prune* is True, the tree is pruned to *maxdepth*, or if that is 0,
        to the value of the *maxdepth* option on the *toctree* node.
        If *titles_only* is True, only toplevel document titles will be in the
        resulting tree.
        If *collapse* is True, all branches not containing docname will
        be collapsed.
        """
        return toctree_adapters._resolve_toctree(
            self,
            docname,
            builder,
            toctree,
            prune=prune,
            maxdepth=maxdepth,
            titles_only=titles_only,
            collapse=collapse,
            includehidden=includehidden,
        )

    def resolve_references(
        self, doctree: nodes.document, fromdocname: str, builder: Builder
    ) -> None:
        self.apply_post_transforms(doctree, fromdocname)

    def apply_post_transforms(self, doctree: nodes.document, docname: str) -> None:
        """Apply all post-transforms."""
        try:
            # set env.docname during applying post-transforms
            backup = copy(self.temp_data)
            self.temp_data['docname'] = docname

            transformer = SphinxTransformer(doctree)
            transformer.set_environment(self)
            transformer.add_transforms(self.app.registry.get_post_transforms())
            transformer.apply_transforms()
        finally:
            self.temp_data = backup

        # allow custom references to be resolved
        self.events.emit('doctree-resolved', doctree, docname)

    def collect_relations(self) -> dict[str, list[str | None]]:
        traversed: set[str] = set()

        relations = {}
        docnames = _traverse_toctree(
            traversed, None, self.config.root_doc, self.toctree_includes
        )
        prev_doc = None
        parent, docname = next(docnames)
        for next_parent, next_doc in docnames:
            relations[docname] = [parent, prev_doc, next_doc]
            prev_doc = docname
            docname = next_doc
            parent = next_parent

        relations[docname] = [parent, prev_doc, None]

        return relations

    def check_consistency(self) -> None:
        """Do consistency checks."""
        included = set().union(*self.included.values())
        for docname in sorted(self.all_docs):
            if docname not in self.files_to_rebuild:
                if docname == self.config.root_doc:
                    # the master file is not included anywhere ;)
                    continue
                if docname in included:
                    # the document is included from other documents
                    continue
                if 'orphan' in self.metadata[docname]:
                    continue
                logger.warning(
                    __("document isn't included in any toctree"), location=docname
                )
        # Call _check_toc_parents here rather than in  _get_toctree_ancestors()
        # because that method is called multiple times per document and would
        # lead to duplicate warnings.
        _check_toc_parents(self.toctree_includes)

        # call check-consistency for all extensions
        self.domains._check_consistency()
        self.events.emit('env-check-consistency', self)


def _differing_config_keys(old: Config, new: Config) -> frozenset[str]:
    """Return a set of keys that differ between two config objects."""
    old_vals = {c.name: c.value for c in old}
    new_vals = {c.name: c.value for c in new}
    not_in_both = old_vals.keys() ^ new_vals.keys()
    different_values = {
        key
        for key in old_vals.keys() & new_vals.keys()
        if stable_str(old_vals[key]) != stable_str(new_vals[key])
    }
    return frozenset(not_in_both | different_values)


def _traverse_toctree(
    traversed: set[str],
    parent: str | None,
    docname: str,
    toctree_includes: dict[str, list[str]],
) -> Iterator[tuple[str | None, str]]:
    if parent == docname:
        logger.warning(
            __('self referenced toctree found. Ignored.'),
            location=docname,
            type='toc',
            subtype='circular',
        )
        return

    # traverse toctree by pre-order
    yield parent, docname
    traversed.add(docname)

    for child in toctree_includes.get(docname, ()):
        for sub_parent, sub_docname in _traverse_toctree(
            traversed, docname, child, toctree_includes
        ):
            if sub_docname not in traversed:
                yield sub_parent, sub_docname
                traversed.add(sub_docname)


def _check_toc_parents(toctree_includes: dict[str, list[str]]) -> None:
    toc_parents: dict[str, list[str]] = {}
    for parent, children in toctree_includes.items():
        for child in children:
            toc_parents.setdefault(child, []).append(parent)

    for doc, parents in sorted(toc_parents.items()):
        if len(parents) > 1:
            logger.info(
                __(
                    'document is referenced in multiple toctrees: %s, selecting: %s <- %s'
                ),
                parents,
                max(parents),
                doc,
                location=doc,
                type='toc',
                subtype='multiple_toc_parents',
            )
