"""
    sphinx.builders.html
    ~~~~~~~~~~~~~~~~~~~~

    Several HTML builders.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import html
import os
import posixpath
import re
import sys
from datetime import datetime
from os import path
from typing import IO, Any, Dict, Iterable, Iterator, List, Set, Tuple, Type
from urllib.parse import quote

from docutils import nodes
from docutils.core import publish_parts
from docutils.frontend import OptionParser
from docutils.io import DocTreeInput, StringOutput
from docutils.nodes import Node
from docutils.utils import relative_path

from sphinx import __display_version__, package_dir
from sphinx import version_info as sphinx_version
from sphinx.application import Sphinx
from sphinx.builders import Builder
from sphinx.config import ENUM, Config
from sphinx.domains import Domain, Index, IndexEntry
from sphinx.environment.adapters.asset import ImageAdapter
from sphinx.environment.adapters.indexentries import IndexEntries
from sphinx.environment.adapters.toctree import TocTree
from sphinx.errors import ConfigError, ThemeError
from sphinx.highlighting import PygmentsBridge
from sphinx.locale import _, __
from sphinx.search import js_index
from sphinx.theming import HTMLThemeFactory
from sphinx.util import isurl, logging, md5, progress_message, status_iterator
from sphinx.util.docutils import is_html5_writer_available, new_document
from sphinx.util.fileutil import copy_asset
from sphinx.util.i18n import format_date
from sphinx.util.inventory import InventoryFile
from sphinx.util.matching import DOTFILES, Matcher, patmatch
from sphinx.util.osutil import copyfile, ensuredir, os_path, relative_uri
from sphinx.util.tags import Tags
from sphinx.writers.html import HTMLTranslator, HTMLWriter

# HTML5 Writer is available or not
if is_html5_writer_available():
    from sphinx.writers.html5 import HTML5Translator
    html5_ready = True
else:
    html5_ready = False

#: the filename for the inventory of objects
INVENTORY_FILENAME = 'objects.inv'

logger = logging.getLogger(__name__)
return_codes_re = re.compile('[\r\n]+')


def get_stable_hash(obj: Any) -> str:
    """
    Return a stable hash for a Python data structure.  We can't just use
    the md5 of str(obj) since for example dictionary items are enumerated
    in unpredictable order due to hash randomization in newer Pythons.
    """
    if isinstance(obj, dict):
        return get_stable_hash(list(obj.items()))
    elif isinstance(obj, (list, tuple)):
        obj = sorted(get_stable_hash(o) for o in obj)
    return md5(str(obj).encode()).hexdigest()


class Stylesheet(str):
    """A metadata of stylesheet.

    To keep compatibility with old themes, an instance of stylesheet behaves as
    its filename (str).
    """

    attributes: Dict[str, str] = None
    filename: str = None
    priority: int = None

    def __new__(cls, filename: str, *args: str, priority: int = 500, **attributes: Any
                ) -> "Stylesheet":
        self = str.__new__(cls, filename)
        self.filename = filename
        self.priority = priority
        self.attributes = attributes
        self.attributes.setdefault('rel', 'stylesheet')
        self.attributes.setdefault('type', 'text/css')
        if args:  # old style arguments (rel, title)
            self.attributes['rel'] = args[0]
            self.attributes['title'] = args[1]

        return self


class JavaScript(str):
    """A metadata of javascript file.

    To keep compatibility with old themes, an instance of javascript behaves as
    its filename (str).
    """

    attributes: Dict[str, str] = None
    filename: str = None
    priority: int = None

    def __new__(cls, filename: str, priority: int = 500, **attributes: str) -> "JavaScript":
        self = str.__new__(cls, filename)
        self.filename = filename
        self.priority = priority
        self.attributes = attributes

        return self


class BuildInfo:
    """buildinfo file manipulator.

    HTMLBuilder and its family are storing their own envdata to ``.buildinfo``.
    This class is a manipulator for the file.
    """

    @classmethod
    def load(cls, f: IO) -> "BuildInfo":
        try:
            lines = f.readlines()
            assert lines[0].rstrip() == '# Sphinx build info version 1'
            assert lines[2].startswith('config: ')
            assert lines[3].startswith('tags: ')

            build_info = BuildInfo()
            build_info.config_hash = lines[2].split()[1].strip()
            build_info.tags_hash = lines[3].split()[1].strip()
            return build_info
        except Exception as exc:
            raise ValueError(__('build info file is broken: %r') % exc) from exc

    def __init__(self, config: Config = None, tags: Tags = None, config_categories: List[str] = []) -> None:  # NOQA
        self.config_hash = ''
        self.tags_hash = ''

        if config:
            values = {c.name: c.value for c in config.filter(config_categories)}
            self.config_hash = get_stable_hash(values)

        if tags:
            self.tags_hash = get_stable_hash(sorted(tags))

    def __eq__(self, other: "BuildInfo") -> bool:  # type: ignore
        return (self.config_hash == other.config_hash and
                self.tags_hash == other.tags_hash)

    def dump(self, f: IO) -> None:
        f.write('# Sphinx build info version 1\n'
                '# This file hashes the configuration used when building these files.'
                ' When it is not found, a full rebuild will be done.\n'
                'config: %s\n'
                'tags: %s\n' %
                (self.config_hash, self.tags_hash))


class StandaloneHTMLBuilder(Builder):
    """
    Builds standalone HTML docs.
    """
    name = 'html'
    format = 'html'
    epilog = __('The HTML pages are in %(outdir)s.')

    copysource = True
    allow_parallel = True
    out_suffix = '.html'
    link_suffix = '.html'  # defaults to matching out_suffix
    indexer_format: Any = js_index
    indexer_dumps_unicode = True
    # create links to original images from images [True/False]
    html_scaled_image_link = True
    supported_image_types = ['image/svg+xml', 'image/png',
                             'image/gif', 'image/jpeg']
    supported_remote_images = True
    supported_data_uri_images = True
    searchindex_filename = 'searchindex.js'
    add_permalinks = True
    allow_sharp_as_current_path = True
    embedded = False  # for things like HTML help or Qt help: suppresses sidebar
    search = True  # for things like HTML help and Apple help: suppress search
    use_index = False
    download_support = True  # enable download role

    imgpath: str = None
    domain_indices: List[Tuple[str, Type[Index], List[Tuple[str, List[IndexEntry]]], bool]] = []  # NOQA

    def __init__(self, app: Sphinx) -> None:
        super().__init__(app)

        # CSS files
        self.css_files: List[Stylesheet] = []

        # JS files
        self.script_files: List[JavaScript] = []

    def init(self) -> None:
        self.build_info = self.create_build_info()
        # basename of images directory
        self.imagedir = '_images'
        # section numbers for headings in the currently visited document
        self.secnumbers: Dict[str, Tuple[int, ...]] = {}
        # currently written docname
        self.current_docname: str = None

        self.init_templates()
        self.init_highlighter()
        self.init_css_files()
        self.init_js_files()

        html_file_suffix = self.get_builder_config('file_suffix', 'html')
        if html_file_suffix is not None:
            self.out_suffix = html_file_suffix

        html_link_suffix = self.get_builder_config('link_suffix', 'html')
        if html_link_suffix is not None:
            self.link_suffix = html_link_suffix
        else:
            self.link_suffix = self.out_suffix

        self.use_index = self.get_builder_config('use_index', 'html')

    def create_build_info(self) -> BuildInfo:
        return BuildInfo(self.config, self.tags, ['html'])

    def _get_translations_js(self) -> str:
        candidates = [path.join(dir, self.config.language,
                                'LC_MESSAGES', 'sphinx.js')
                      for dir in self.config.locale_dirs] + \
                     [path.join(package_dir, 'locale', self.config.language,
                                'LC_MESSAGES', 'sphinx.js'),
                      path.join(sys.prefix, 'share/sphinx/locale',
                                self.config.language, 'sphinx.js')]

        for jsfile in candidates:
            if path.isfile(jsfile):
                return jsfile
        return None

    def _get_style_filename(self) -> str:
        if self.config.html_style is not None:
            return self.config.html_style
        elif self.theme:
            return self.theme.get_config('theme', 'stylesheet')
        else:
            return 'default.css'

    def get_theme_config(self) -> Tuple[str, Dict]:
        return self.config.html_theme, self.config.html_theme_options

    def init_templates(self) -> None:
        theme_factory = HTMLThemeFactory(self.app)
        themename, themeoptions = self.get_theme_config()
        self.theme = theme_factory.create(themename)
        self.theme_options = themeoptions.copy()
        self.create_template_bridge()
        self.templates.init(self, self.theme)

    def init_highlighter(self) -> None:
        # determine Pygments style and create the highlighter
        if self.config.pygments_style is not None:
            style = self.config.pygments_style
        elif self.theme:
            style = self.theme.get_config('theme', 'pygments_style', 'none')
        else:
            style = 'sphinx'
        self.highlighter = PygmentsBridge('html', style)

        if self.theme:
            dark_style = self.theme.get_config('theme', 'pygments_dark_style', None)
        else:
            dark_style = None

        if dark_style is not None:
            self.dark_highlighter = PygmentsBridge('html', dark_style)
            self.app.add_css_file('pygments_dark.css',
                                  media='(prefers-color-scheme: dark)',
                                  id='pygments_dark_css')
        else:
            self.dark_highlighter = None

    def init_css_files(self) -> None:
        self.css_files = []
        self.add_css_file('pygments.css', priority=200)
        self.add_css_file(self._get_style_filename(), priority=200)

        for filename, attrs in self.app.registry.css_files:
            self.add_css_file(filename, **attrs)

        for filename, attrs in self.get_builder_config('css_files', 'html'):
            attrs.setdefault('priority', 800)  # User's CSSs are loaded after extensions'
            self.add_css_file(filename, **attrs)

    def add_css_file(self, filename: str, **kwargs: Any) -> None:
        if '://' not in filename:
            filename = posixpath.join('_static', filename)

        self.css_files.append(Stylesheet(filename, **kwargs))

    def init_js_files(self) -> None:
        self.script_files = []
        self.add_js_file('documentation_options.js', id="documentation_options",
                         data_url_root='', priority=200)
        self.add_js_file('jquery.js', priority=200)
        self.add_js_file('underscore.js', priority=200)
        self.add_js_file('doctools.js', priority=200)

        for filename, attrs in self.app.registry.js_files:
            self.add_js_file(filename, **attrs)

        for filename, attrs in self.get_builder_config('js_files', 'html'):
            attrs.setdefault('priority', 800)  # User's JSs are loaded after extensions'
            self.add_js_file(filename, **attrs)

        if self.config.language and self._get_translations_js():
            self.add_js_file('translations.js')

    def add_js_file(self, filename: str, **kwargs: Any) -> None:
        if filename and '://' not in filename:
            filename = posixpath.join('_static', filename)

        self.script_files.append(JavaScript(filename, **kwargs))

    @property
    def default_translator_class(self) -> Type[nodes.NodeVisitor]:  # type: ignore
        if not html5_ready or self.config.html4_writer:
            return HTMLTranslator
        else:
            return HTML5Translator

    @property
    def math_renderer_name(self) -> str:
        name = self.get_builder_config('math_renderer', 'html')
        if name is not None:
            # use given name
            return name
        else:
            # not given: choose a math_renderer from registered ones as possible
            renderers = list(self.app.registry.html_inline_math_renderers)
            if len(renderers) == 1:
                # only default math_renderer (mathjax) is registered
                return renderers[0]
            elif len(renderers) == 2:
                # default and another math_renderer are registered; prior the another
                renderers.remove('mathjax')
                return renderers[0]
            else:
                # many math_renderers are registered. can't choose automatically!
                return None

    def get_outdated_docs(self) -> Iterator[str]:
        try:
            with open(path.join(self.outdir, '.buildinfo')) as fp:
                buildinfo = BuildInfo.load(fp)

            if self.build_info != buildinfo:
                logger.debug('[build target] did not match: build_info ')
                yield from self.env.found_docs
                return
        except ValueError as exc:
            logger.warning(__('Failed to read build info file: %r'), exc)
        except OSError:
            # ignore errors on reading
            pass

        if self.templates:
            template_mtime = self.templates.newest_template_mtime()
        else:
            template_mtime = 0
        for docname in self.env.found_docs:
            if docname not in self.env.all_docs:
                logger.debug('[build target] did not in env: %r', docname)
                yield docname
                continue
            targetname = self.get_outfilename(docname)
            try:
                targetmtime = path.getmtime(targetname)
            except Exception:
                targetmtime = 0
            try:
                srcmtime = max(path.getmtime(self.env.doc2path(docname)),
                               template_mtime)
                if srcmtime > targetmtime:
                    logger.debug(
                        '[build target] targetname %r(%s), template(%s), docname %r(%s)',
                        targetname,
                        datetime.utcfromtimestamp(targetmtime),
                        datetime.utcfromtimestamp(template_mtime),
                        docname,
                        datetime.utcfromtimestamp(path.getmtime(self.env.doc2path(docname))),
                    )
                    yield docname
            except OSError:
                # source doesn't exist anymore
                pass

    def get_asset_paths(self) -> List[str]:
        return self.config.html_extra_path + self.config.html_static_path

    def render_partial(self, node: Node) -> Dict[str, str]:
        """Utility: Render a lone doctree node."""
        if node is None:
            return {'fragment': ''}
        doc = new_document('<partial node>')
        doc.append(node)

        writer = HTMLWriter(self)
        return publish_parts(reader_name='doctree',
                             writer=writer,
                             source_class=DocTreeInput,
                             settings_overrides={'output_encoding': 'unicode'},
                             source=doc)

    def prepare_writing(self, docnames: Set[str]) -> None:
        # create the search indexer
        self.indexer = None
        if self.search:
            from sphinx.search import IndexBuilder
            lang = self.config.html_search_language or self.config.language
            if not lang:
                lang = 'en'
            self.indexer = IndexBuilder(self.env, lang,
                                        self.config.html_search_options,
                                        self.config.html_search_scorer)
            self.load_indexer(docnames)

        self.docwriter = HTMLWriter(self)
        self.docsettings: Any = OptionParser(
            defaults=self.env.settings,
            components=(self.docwriter,),
            read_config_files=True).get_default_values()
        self.docsettings.compact_lists = bool(self.config.html_compact_lists)

        # determine the additional indices to include
        self.domain_indices = []
        # html_domain_indices can be False/True or a list of index names
        indices_config = self.config.html_domain_indices
        if indices_config:
            for domain_name in sorted(self.env.domains):
                domain: Domain = self.env.domains[domain_name]
                for indexcls in domain.indices:
                    indexname = '%s-%s' % (domain.name, indexcls.name)
                    if isinstance(indices_config, list):
                        if indexname not in indices_config:
                            continue
                    content, collapse = indexcls(domain).generate()
                    if content:
                        self.domain_indices.append(
                            (indexname, indexcls, content, collapse))

        # format the "last updated on" string, only once is enough since it
        # typically doesn't include the time of day
        lufmt = self.config.html_last_updated_fmt
        if lufmt is not None:
            self.last_updated = format_date(lufmt or _('%b %d, %Y'),
                                            language=self.config.language)
        else:
            self.last_updated = None

        # If the logo or favicon are urls, keep them as-is, otherwise
        # strip the relative path as the files will be copied into _static.
        logo = self.config.html_logo or ''
        favicon = self.config.html_favicon or ''

        if not isurl(logo):
            logo = path.basename(logo)
        if not isurl(favicon):
            favicon = path.basename(favicon)

        self.relations = self.env.collect_relations()

        rellinks: List[Tuple[str, str, str, str]] = []
        if self.use_index:
            rellinks.append(('genindex', _('General Index'), 'I', _('index')))
        for indexname, indexcls, _content, _collapse in self.domain_indices:
            # if it has a short name
            if indexcls.shortname:
                rellinks.append((indexname, indexcls.localname,
                                 '', indexcls.shortname))

        # back up script_files and css_files to allow adding JS/CSS files to a specific page.
        self._script_files = list(self.script_files)
        self._css_files = list(self.css_files)

        self.globalcontext = {
            'embedded': self.embedded,
            'project': self.config.project,
            'release': return_codes_re.sub('', self.config.release),
            'version': self.config.version,
            'last_updated': self.last_updated,
            'copyright': self.config.copyright,
            'master_doc': self.config.root_doc,
            'root_doc': self.config.root_doc,
            'use_opensearch': self.config.html_use_opensearch,
            'docstitle': self.config.html_title,
            'shorttitle': self.config.html_short_title,
            'show_copyright': self.config.html_show_copyright,
            'show_sphinx': self.config.html_show_sphinx,
            'has_source': self.config.html_copy_source,
            'show_source': self.config.html_show_sourcelink,
            'sourcelink_suffix': self.config.html_sourcelink_suffix,
            'file_suffix': self.out_suffix,
            'link_suffix': self.link_suffix,
            'script_files': self.script_files,
            'language': self.config.language,
            'css_files': self.css_files,
            'sphinx_version': __display_version__,
            'sphinx_version_tuple': sphinx_version,
            'style': self._get_style_filename(),
            'rellinks': rellinks,
            'builder': self.name,
            'parents': [],
            'logo': logo,
            'favicon': favicon,
            'html5_doctype': html5_ready and not self.config.html4_writer,
        }
        if self.theme:
            self.globalcontext.update(
                ('theme_' + key, val) for (key, val) in
                self.theme.get_options(self.theme_options).items())
        self.globalcontext.update(self.config.html_context)

    def get_doc_context(self, docname: str, body: str, metatags: str) -> Dict[str, Any]:
        """Collect items for the template context of a page."""
        # find out relations
        prev = next = None
        parents = []
        rellinks = self.globalcontext['rellinks'][:]
        related = self.relations.get(docname)
        titles = self.env.titles
        if related and related[2]:
            try:
                next = {
                    'link': self.get_relative_uri(docname, related[2]),
                    'title': self.render_partial(titles[related[2]])['title']
                }
                rellinks.append((related[2], next['title'], 'N', _('next')))
            except KeyError:
                next = None
        if related and related[1]:
            try:
                prev = {
                    'link': self.get_relative_uri(docname, related[1]),
                    'title': self.render_partial(titles[related[1]])['title']
                }
                rellinks.append((related[1], prev['title'], 'P', _('previous')))
            except KeyError:
                # the relation is (somehow) not in the TOC tree, handle
                # that gracefully
                prev = None
        while related and related[0]:
            try:
                parents.append(
                    {'link': self.get_relative_uri(docname, related[0]),
                     'title': self.render_partial(titles[related[0]])['title']})
            except KeyError:
                pass
            related = self.relations.get(related[0])
        if parents:
            # remove link to the master file; we have a generic
            # "back to index" link already
            parents.pop()
        parents.reverse()

        # title rendered as HTML
        title_node = self.env.longtitles.get(docname)
        title = self.render_partial(title_node)['title'] if title_node else ''

        # Suffix for the document
        source_suffix = self.env.doc2path(docname, False)[len(docname):]

        # the name for the copied source
        if self.config.html_copy_source:
            sourcename = docname + source_suffix
            if source_suffix != self.config.html_sourcelink_suffix:
                sourcename += self.config.html_sourcelink_suffix
        else:
            sourcename = ''

        # metadata for the document
        meta = self.env.metadata.get(docname)

        # local TOC and global TOC tree
        self_toc = TocTree(self.env).get_toc_for(docname, self)
        toc = self.render_partial(self_toc)['fragment']

        return {
            'parents': parents,
            'prev': prev,
            'next': next,
            'title': title,
            'meta': meta,
            'body': body,
            'metatags': metatags,
            'rellinks': rellinks,
            'sourcename': sourcename,
            'toc': toc,
            # only display a TOC if there's more than one item to show
            'display_toc': (self.env.toc_num_entries[docname] > 1),
            'page_source_suffix': source_suffix,
        }

    def write_doc(self, docname: str, doctree: nodes.document) -> None:
        destination = StringOutput(encoding='utf-8')
        doctree.settings = self.docsettings

        self.secnumbers = self.env.toc_secnumbers.get(docname, {})
        self.fignumbers = self.env.toc_fignumbers.get(docname, {})
        self.imgpath = relative_uri(self.get_target_uri(docname), '_images')
        self.dlpath = relative_uri(self.get_target_uri(docname), '_downloads')
        self.current_docname = docname
        self.docwriter.write(doctree, destination)
        self.docwriter.assemble_parts()
        body = self.docwriter.parts['fragment']
        metatags = self.docwriter.clean_meta

        ctx = self.get_doc_context(docname, body, metatags)
        self.handle_page(docname, ctx, event_arg=doctree)

    def write_doc_serialized(self, docname: str, doctree: nodes.document) -> None:
        self.imgpath = relative_uri(self.get_target_uri(docname), self.imagedir)
        self.post_process_images(doctree)
        title_node = self.env.longtitles.get(docname)
        title = self.render_partial(title_node)['title'] if title_node else ''
        self.index_page(docname, doctree, title)

    def finish(self) -> None:
        self.finish_tasks.add_task(self.gen_indices)
        self.finish_tasks.add_task(self.gen_pages_from_extensions)
        self.finish_tasks.add_task(self.gen_additional_pages)
        self.finish_tasks.add_task(self.copy_image_files)
        self.finish_tasks.add_task(self.copy_download_files)
        self.finish_tasks.add_task(self.copy_static_files)
        self.finish_tasks.add_task(self.copy_extra_files)
        self.finish_tasks.add_task(self.write_buildinfo)

        # dump the search index
        self.handle_finish()

    @progress_message(__('generating indices'))
    def gen_indices(self) -> None:
        # the global general index
        if self.use_index:
            self.write_genindex()

        # the global domain-specific indices
        self.write_domain_indices()

    def gen_pages_from_extensions(self) -> None:
        # pages from extensions
        for pagelist in self.events.emit('html-collect-pages'):
            for pagename, context, template in pagelist:
                self.handle_page(pagename, context, template)

    @progress_message(__('writing additional pages'))
    def gen_additional_pages(self) -> None:
        # additional pages from conf.py
        for pagename, template in self.config.html_additional_pages.items():
            logger.info(pagename + ' ', nonl=True)
            self.handle_page(pagename, {}, template)

        # the search page
        if self.search:
            logger.info('search ', nonl=True)
            self.handle_page('search', {}, 'search.html')

        # the opensearch xml file
        if self.config.html_use_opensearch and self.search:
            logger.info('opensearch ', nonl=True)
            fn = path.join(self.outdir, '_static', 'opensearch.xml')
            self.handle_page('opensearch', {}, 'opensearch.xml', outfilename=fn)

    def write_genindex(self) -> None:
        # the total count of lines for each index letter, used to distribute
        # the entries into two columns
        genindex = IndexEntries(self.env).create_index(self)
        indexcounts = []
        for _k, entries in genindex:
            indexcounts.append(sum(1 + len(subitems)
                                   for _, (_, subitems, _) in entries))

        genindexcontext = {
            'genindexentries': genindex,
            'genindexcounts': indexcounts,
            'split_index': self.config.html_split_index,
        }
        logger.info('genindex ', nonl=True)

        if self.config.html_split_index:
            self.handle_page('genindex', genindexcontext,
                             'genindex-split.html')
            self.handle_page('genindex-all', genindexcontext,
                             'genindex.html')
            for (key, entries), count in zip(genindex, indexcounts):
                ctx = {'key': key, 'entries': entries, 'count': count,
                       'genindexentries': genindex}
                self.handle_page('genindex-' + key, ctx,
                                 'genindex-single.html')
        else:
            self.handle_page('genindex', genindexcontext, 'genindex.html')

    def write_domain_indices(self) -> None:
        for indexname, indexcls, content, collapse in self.domain_indices:
            indexcontext = {
                'indextitle': indexcls.localname,
                'content': content,
                'collapse_index': collapse,
            }
            logger.info(indexname + ' ', nonl=True)
            self.handle_page(indexname, indexcontext, 'domainindex.html')

    def copy_image_files(self) -> None:
        if self.images:
            stringify_func = ImageAdapter(self.app.env).get_original_image_uri
            ensuredir(path.join(self.outdir, self.imagedir))
            for src in status_iterator(self.images, __('copying images... '), "brown",
                                       len(self.images), self.app.verbosity,
                                       stringify_func=stringify_func):
                dest = self.images[src]
                try:
                    copyfile(path.join(self.srcdir, src),
                             path.join(self.outdir, self.imagedir, dest))
                except Exception as err:
                    logger.warning(__('cannot copy image file %r: %s'),
                                   path.join(self.srcdir, src), err)

    def copy_download_files(self) -> None:
        def to_relpath(f: str) -> str:
            return relative_path(self.srcdir, f)

        # copy downloadable files
        if self.env.dlfiles:
            ensuredir(path.join(self.outdir, '_downloads'))
            for src in status_iterator(self.env.dlfiles, __('copying downloadable files... '),
                                       "brown", len(self.env.dlfiles), self.app.verbosity,
                                       stringify_func=to_relpath):
                try:
                    dest = path.join(self.outdir, '_downloads', self.env.dlfiles[src][1])
                    ensuredir(path.dirname(dest))
                    copyfile(path.join(self.srcdir, src), dest)
                except OSError as err:
                    logger.warning(__('cannot copy downloadable file %r: %s'),
                                   path.join(self.srcdir, src), err)

    def create_pygments_style_file(self) -> None:
        """create a style file for pygments."""
        with open(path.join(self.outdir, '_static', 'pygments.css'), 'w') as f:
            f.write(self.highlighter.get_stylesheet())

        if self.dark_highlighter:
            with open(path.join(self.outdir, '_static', 'pygments_dark.css'), 'w') as f:
                f.write(self.dark_highlighter.get_stylesheet())

    def copy_translation_js(self) -> None:
        """Copy a JavaScript file for translations."""
        if self.config.language is not None:
            jsfile = self._get_translations_js()
            if jsfile:
                copyfile(jsfile, path.join(self.outdir, '_static', 'translations.js'))

    def copy_stemmer_js(self) -> None:
        """Copy a JavaScript file for stemmer."""
        if self.indexer is not None:
            if hasattr(self.indexer, 'get_js_stemmer_rawcodes'):
                for jsfile in self.indexer.get_js_stemmer_rawcodes():
                    copyfile(jsfile, path.join(self.outdir, '_static', path.basename(jsfile)))
            else:
                jsfile = self.indexer.get_js_stemmer_rawcode()
                if jsfile:
                    copyfile(jsfile, path.join(self.outdir, '_static', '_stemmer.js'))

    def copy_theme_static_files(self, context: Dict) -> None:
        def onerror(filename: str, error: Exception) -> None:
            logger.warning(__('Failed to copy a file in html_static_file: %s: %r'),
                           filename, error)

        if self.theme:
            for entry in self.theme.get_theme_dirs()[::-1]:
                copy_asset(path.join(entry, 'static'),
                           path.join(self.outdir, '_static'),
                           excluded=DOTFILES, context=context,
                           renderer=self.templates, onerror=onerror)

    def copy_html_static_files(self, context: Dict) -> None:
        def onerror(filename: str, error: Exception) -> None:
            logger.warning(__('Failed to copy a file in html_static_file: %s: %r'),
                           filename, error)

        excluded = Matcher(self.config.exclude_patterns + ["**/.*"])
        for entry in self.config.html_static_path:
            copy_asset(path.join(self.confdir, entry),
                       path.join(self.outdir, '_static'),
                       excluded, context=context, renderer=self.templates, onerror=onerror)

    def copy_html_logo(self) -> None:
        if self.config.html_logo and not isurl(self.config.html_logo):
            copy_asset(path.join(self.confdir, self.config.html_logo),
                       path.join(self.outdir, '_static'))

    def copy_html_favicon(self) -> None:
        if self.config.html_favicon and not isurl(self.config.html_favicon):
            copy_asset(path.join(self.confdir, self.config.html_favicon),
                       path.join(self.outdir, '_static'))

    def copy_static_files(self) -> None:
        try:
            with progress_message(__('copying static files')):
                ensuredir(path.join(self.outdir, '_static'))

                # prepare context for templates
                context = self.globalcontext.copy()
                if self.indexer is not None:
                    context.update(self.indexer.context_for_searchtool())

                self.create_pygments_style_file()
                self.copy_translation_js()
                self.copy_stemmer_js()
                self.copy_theme_static_files(context)
                self.copy_html_static_files(context)
                self.copy_html_logo()
                self.copy_html_favicon()
        except OSError as err:
            logger.warning(__('cannot copy static file %r'), err)

    def copy_extra_files(self) -> None:
        """copy html_extra_path files."""
        try:
            with progress_message(__('copying extra files')):
                excluded = Matcher(self.config.exclude_patterns)
                for extra_path in self.config.html_extra_path:
                    entry = path.join(self.confdir, extra_path)
                    copy_asset(entry, self.outdir, excluded)
        except OSError as err:
            logger.warning(__('cannot copy extra file %r'), err)

    def write_buildinfo(self) -> None:
        try:
            with open(path.join(self.outdir, '.buildinfo'), 'w') as fp:
                self.build_info.dump(fp)
        except OSError as exc:
            logger.warning(__('Failed to write build info file: %r'), exc)

    def cleanup(self) -> None:
        # clean up theme stuff
        if self.theme:
            self.theme.cleanup()

    def post_process_images(self, doctree: Node) -> None:
        """Pick the best candidate for an image and link down-scaled images to
        their high res version.
        """
        Builder.post_process_images(self, doctree)

        if self.config.html_scaled_image_link and self.html_scaled_image_link:
            for node in doctree.findall(nodes.image):
                if not any((key in node) for key in ['scale', 'width', 'height']):
                    # resizing options are not given. scaled image link is available
                    # only for resized images.
                    continue
                elif isinstance(node.parent, nodes.reference):
                    # A image having hyperlink target
                    continue
                elif 'no-scaled-link' in node['classes']:
                    # scaled image link is disabled for this node
                    continue

                uri = node['uri']
                reference = nodes.reference('', '', internal=True)
                if uri in self.images:
                    reference['refuri'] = posixpath.join(self.imgpath,
                                                         self.images[uri])
                else:
                    reference['refuri'] = uri
                node.replace_self(reference)
                reference.append(node)

    def load_indexer(self, docnames: Iterable[str]) -> None:
        keep = set(self.env.all_docs) - set(docnames)
        try:
            searchindexfn = path.join(self.outdir, self.searchindex_filename)
            if self.indexer_dumps_unicode:
                with open(searchindexfn, encoding='utf-8') as ft:
                    self.indexer.load(ft, self.indexer_format)
            else:
                with open(searchindexfn, 'rb') as fb:
                    self.indexer.load(fb, self.indexer_format)
        except (OSError, ValueError):
            if keep:
                logger.warning(__('search index couldn\'t be loaded, but not all '
                                  'documents will be built: the index will be '
                                  'incomplete.'))
        # delete all entries for files that will be rebuilt
        self.indexer.prune(keep)

    def index_page(self, pagename: str, doctree: nodes.document, title: str) -> None:
        # only index pages with title
        if self.indexer is not None and title:
            filename = self.env.doc2path(pagename, base=None)
            metadata = self.env.metadata.get(pagename, {})
            if 'nosearch' in metadata:
                self.indexer.feed(pagename, filename, '', new_document(''))
            else:
                self.indexer.feed(pagename, filename, title, doctree)

    def _get_local_toctree(self, docname: str, collapse: bool = True, **kwargs: Any) -> str:
        if 'includehidden' not in kwargs:
            kwargs['includehidden'] = False
        if kwargs.get('maxdepth') == '':
            kwargs.pop('maxdepth')
        return self.render_partial(TocTree(self.env).get_toctree_for(
            docname, self, collapse, **kwargs))['fragment']

    def get_outfilename(self, pagename: str) -> str:
        return path.join(self.outdir, os_path(pagename) + self.out_suffix)

    def add_sidebars(self, pagename: str, ctx: Dict) -> None:
        def has_wildcard(pattern: str) -> bool:
            return any(char in pattern for char in '*?[')

        sidebars = None
        matched = None
        customsidebar = None

        # default sidebars settings for selected theme
        if self.theme.name == 'alabaster':
            # provide default settings for alabaster (for compatibility)
            # Note: this will be removed before Sphinx-2.0
            try:
                # get default sidebars settings from alabaster (if defined)
                theme_default_sidebars = self.theme.config.get('theme', 'sidebars')
                if theme_default_sidebars:
                    sidebars = [name.strip() for name in theme_default_sidebars.split(',')]
            except Exception:
                # fallback to better default settings
                sidebars = ['about.html', 'navigation.html', 'relations.html',
                            'searchbox.html', 'donate.html']
        else:
            theme_default_sidebars = self.theme.get_config('theme', 'sidebars', None)
            if theme_default_sidebars:
                sidebars = [name.strip() for name in theme_default_sidebars.split(',')]

        # user sidebar settings
        html_sidebars = self.get_builder_config('sidebars', 'html')
        for pattern, patsidebars in html_sidebars.items():
            if patmatch(pagename, pattern):
                if matched:
                    if has_wildcard(pattern):
                        # warn if both patterns contain wildcards
                        if has_wildcard(matched):
                            logger.warning(__('page %s matches two patterns in '
                                              'html_sidebars: %r and %r'),
                                           pagename, matched, pattern)
                        # else the already matched pattern is more specific
                        # than the present one, because it contains no wildcard
                        continue
                matched = pattern
                sidebars = patsidebars

        if sidebars is None:
            # keep defaults
            pass

        ctx['sidebars'] = sidebars
        ctx['customsidebar'] = customsidebar

    # --------- these are overwritten by the serialization builder

    def get_target_uri(self, docname: str, typ: str = None) -> str:
        return quote(docname) + self.link_suffix

    def handle_page(self, pagename: str, addctx: Dict, templatename: str = 'page.html',
                    outfilename: str = None, event_arg: Any = None) -> None:
        ctx = self.globalcontext.copy()
        # current_page_name is backwards compatibility
        ctx['pagename'] = ctx['current_page_name'] = pagename
        ctx['encoding'] = self.config.html_output_encoding
        default_baseuri = self.get_target_uri(pagename)
        # in the singlehtml builder, default_baseuri still contains an #anchor
        # part, which relative_uri doesn't really like...
        default_baseuri = default_baseuri.rsplit('#', 1)[0]

        if self.config.html_baseurl:
            ctx['pageurl'] = posixpath.join(self.config.html_baseurl,
                                            pagename + self.out_suffix)
        else:
            ctx['pageurl'] = None

        def pathto(otheruri: str, resource: bool = False, baseuri: str = default_baseuri) -> str:  # NOQA
            if resource and '://' in otheruri:
                # allow non-local resources given by scheme
                return otheruri
            elif not resource:
                otheruri = self.get_target_uri(otheruri)
            uri = relative_uri(baseuri, otheruri) or '#'
            if uri == '#' and not self.allow_sharp_as_current_path:
                uri = baseuri
            return uri
        ctx['pathto'] = pathto

        def hasdoc(name: str) -> bool:
            if name in self.env.all_docs:
                return True
            elif name == 'search' and self.search:
                return True
            elif name == 'genindex' and self.get_builder_config('use_index', 'html'):
                return True
            return False
        ctx['hasdoc'] = hasdoc

        ctx['toctree'] = lambda **kwargs: self._get_local_toctree(pagename, **kwargs)
        self.add_sidebars(pagename, ctx)
        ctx.update(addctx)

        # revert script_files and css_files
        self.script_files[:] = self._script_files
        self.css_files[:] = self._css_files

        self.update_page_context(pagename, templatename, ctx, event_arg)
        newtmpl = self.app.emit_firstresult('html-page-context', pagename,
                                            templatename, ctx, event_arg)
        if newtmpl:
            templatename = newtmpl

        # sort JS/CSS before rendering HTML
        try:
            # Convert script_files to list to support non-list script_files (refs: #8889)
            ctx['script_files'] = sorted(list(ctx['script_files']), key=lambda js: js.priority)
        except AttributeError:
            # Skip sorting if users modifies script_files directly (maybe via `html_context`).
            # refs: #8885
            #
            # Note: priority sorting feature will not work in this case.
            pass

        try:
            ctx['css_files'] = sorted(list(ctx['css_files']), key=lambda css: css.priority)
        except AttributeError:
            pass

        try:
            output = self.templates.render(templatename, ctx)
        except UnicodeError:
            logger.warning(__("a Unicode error occurred when rendering the page %s. "
                              "Please make sure all config values that contain "
                              "non-ASCII content are Unicode strings."), pagename)
            return
        except Exception as exc:
            raise ThemeError(__("An error happened in rendering the page %s.\nReason: %r") %
                             (pagename, exc)) from exc

        if not outfilename:
            outfilename = self.get_outfilename(pagename)
        # outfilename's path is in general different from self.outdir
        ensuredir(path.dirname(outfilename))
        try:
            with open(outfilename, 'w', encoding=ctx['encoding'],
                      errors='xmlcharrefreplace') as f:
                f.write(output)
        except OSError as err:
            logger.warning(__("error writing file %s: %s"), outfilename, err)
        if self.copysource and ctx.get('sourcename'):
            # copy the source file for the "show source" link
            source_name = path.join(self.outdir, '_sources',
                                    os_path(ctx['sourcename']))
            ensuredir(path.dirname(source_name))
            copyfile(self.env.doc2path(pagename), source_name)

    def update_page_context(self, pagename: str, templatename: str,
                            ctx: Dict, event_arg: Any) -> None:
        pass

    def handle_finish(self) -> None:
        if self.indexer:
            self.finish_tasks.add_task(self.dump_search_index)
        self.finish_tasks.add_task(self.dump_inventory)

    @progress_message(__('dumping object inventory'))
    def dump_inventory(self) -> None:
        InventoryFile.dump(path.join(self.outdir, INVENTORY_FILENAME), self.env, self)

    def dump_search_index(self) -> None:
        with progress_message(__('dumping search index in %s') % self.indexer.label()):
            self.indexer.prune(self.env.all_docs)
            searchindexfn = path.join(self.outdir, self.searchindex_filename)
            # first write to a temporary file, so that if dumping fails,
            # the existing index won't be overwritten
            if self.indexer_dumps_unicode:
                with open(searchindexfn + '.tmp', 'w', encoding='utf-8') as ft:
                    self.indexer.dump(ft, self.indexer_format)
            else:
                with open(searchindexfn + '.tmp', 'wb') as fb:
                    self.indexer.dump(fb, self.indexer_format)
            os.replace(searchindexfn + '.tmp', searchindexfn)


def convert_html_css_files(app: Sphinx, config: Config) -> None:
    """This converts string styled html_css_files to tuple styled one."""
    html_css_files: List[Tuple[str, Dict]] = []
    for entry in config.html_css_files:
        if isinstance(entry, str):
            html_css_files.append((entry, {}))
        else:
            try:
                filename, attrs = entry
                html_css_files.append((filename, attrs))
            except Exception:
                logger.warning(__('invalid css_file: %r, ignored'), entry)
                continue

    config.html_css_files = html_css_files  # type: ignore


def convert_html_js_files(app: Sphinx, config: Config) -> None:
    """This converts string styled html_js_files to tuple styled one."""
    html_js_files: List[Tuple[str, Dict]] = []
    for entry in config.html_js_files:
        if isinstance(entry, str):
            html_js_files.append((entry, {}))
        else:
            try:
                filename, attrs = entry
                html_js_files.append((filename, attrs))
            except Exception:
                logger.warning(__('invalid js_file: %r, ignored'), entry)
                continue

    config.html_js_files = html_js_files  # type: ignore


def setup_css_tag_helper(app: Sphinx, pagename: str, templatename: str,
                         context: Dict, doctree: Node) -> None:
    """Set up css_tag() template helper.

    .. note:: This set up function is added to keep compatibility with webhelper.
    """
    pathto = context.get('pathto')

    def css_tag(css: Stylesheet) -> str:
        attrs = []
        for key in sorted(css.attributes):
            value = css.attributes[key]
            if value is not None:
                attrs.append('%s="%s"' % (key, html.escape(value, True)))
        attrs.append('href="%s"' % pathto(css.filename, resource=True))
        return '<link %s />' % ' '.join(attrs)

    context['css_tag'] = css_tag


def setup_js_tag_helper(app: Sphinx, pagename: str, templatename: str,
                        context: Dict, doctree: Node) -> None:
    """Set up js_tag() template helper.

    .. note:: This set up function is added to keep compatibility with webhelper.
    """
    pathto = context.get('pathto')

    def js_tag(js: JavaScript) -> str:
        attrs = []
        body = ''
        if isinstance(js, JavaScript):
            for key in sorted(js.attributes):
                value = js.attributes[key]
                if value is not None:
                    if key == 'body':
                        body = value
                    elif key == 'data_url_root':
                        attrs.append('data-url_root="%s"' % pathto('', resource=True))
                    else:
                        attrs.append('%s="%s"' % (key, html.escape(value, True)))
            if js.filename:
                attrs.append('src="%s"' % pathto(js.filename, resource=True))
        else:
            # str value (old styled)
            attrs.append('src="%s"' % pathto(js, resource=True))

        if attrs:
            return '<script %s>%s</script>' % (' '.join(attrs), body)
        else:
            return '<script>%s</script>' % body

    context['js_tag'] = js_tag


def setup_resource_paths(app: Sphinx, pagename: str, templatename: str,
                         context: Dict, doctree: Node) -> None:
    """Set up relative resource paths."""
    pathto = context.get('pathto')

    # favicon_url
    favicon = context.get('favicon')
    if favicon and not isurl(favicon):
        context['favicon_url'] = pathto('_static/' + favicon, resource=True)
    else:
        context['favicon_url'] = favicon

    # logo_url
    logo = context.get('logo')
    if logo and not isurl(logo):
        context['logo_url'] = pathto('_static/' + logo, resource=True)
    else:
        context['logo_url'] = logo


def validate_math_renderer(app: Sphinx) -> None:
    if app.builder.format != 'html':
        return

    name = app.builder.math_renderer_name  # type: ignore
    if name is None:
        raise ConfigError(__('Many math_renderers are registered. '
                             'But no math_renderer is selected.'))
    elif name not in app.registry.html_inline_math_renderers:
        raise ConfigError(__('Unknown math_renderer %r is given.') % name)


def validate_html_extra_path(app: Sphinx, config: Config) -> None:
    """Check html_extra_paths setting."""
    for entry in config.html_extra_path[:]:
        extra_path = path.normpath(path.join(app.confdir, entry))
        if not path.exists(extra_path):
            logger.warning(__('html_extra_path entry %r does not exist'), entry)
            config.html_extra_path.remove(entry)
        elif (path.splitdrive(app.outdir)[0] == path.splitdrive(extra_path)[0] and
              path.commonpath([app.outdir, extra_path]) == app.outdir):
            logger.warning(__('html_extra_path entry %r is placed inside outdir'), entry)
            config.html_extra_path.remove(entry)


def validate_html_static_path(app: Sphinx, config: Config) -> None:
    """Check html_static_paths setting."""
    for entry in config.html_static_path[:]:
        static_path = path.normpath(path.join(app.confdir, entry))
        if not path.exists(static_path):
            logger.warning(__('html_static_path entry %r does not exist'), entry)
            config.html_static_path.remove(entry)
        elif (path.splitdrive(app.outdir)[0] == path.splitdrive(static_path)[0] and
              path.commonpath([app.outdir, static_path]) == app.outdir):
            logger.warning(__('html_static_path entry %r is placed inside outdir'), entry)
            config.html_static_path.remove(entry)


def validate_html_logo(app: Sphinx, config: Config) -> None:
    """Check html_logo setting."""
    if (config.html_logo and
            not path.isfile(path.join(app.confdir, config.html_logo)) and
            not isurl(config.html_logo)):
        logger.warning(__('logo file %r does not exist'), config.html_logo)
        config.html_logo = None  # type: ignore


def validate_html_favicon(app: Sphinx, config: Config) -> None:
    """Check html_favicon setting."""
    if (config.html_favicon and
            not path.isfile(path.join(app.confdir, config.html_favicon)) and
            not isurl(config.html_favicon)):
        logger.warning(__('favicon file %r does not exist'), config.html_favicon)
        config.html_favicon = None  # type: ignore


class _stable_repr_object():

    def __repr__(self):
        return '<object>'


UNSET = _stable_repr_object()


def migrate_html_add_permalinks(app: Sphinx, config: Config) -> None:
    """Migrate html_add_permalinks to html_permalinks*."""
    html_add_permalinks = config.html_add_permalinks
    if html_add_permalinks is UNSET:
        return

    # RemovedInSphinx60Warning
    logger.warning(__('html_add_permalinks has been deprecated since v3.5.0. '
                      'Please use html_permalinks and html_permalinks_icon instead.'))
    if not html_add_permalinks:
        config.html_permalinks = False  # type: ignore[attr-defined]
        return

    config.html_permalinks_icon = html.escape(  # type: ignore[attr-defined]
        html_add_permalinks
    )

# for compatibility
import sphinxcontrib.serializinghtml  # NOQA

import sphinx.builders.dirhtml  # NOQA
import sphinx.builders.singlehtml  # NOQA


def setup(app: Sphinx) -> Dict[str, Any]:
    # builders
    app.add_builder(StandaloneHTMLBuilder)

    # config values
    app.add_config_value('html_theme', 'alabaster', 'html')
    app.add_config_value('html_theme_path', [], 'html')
    app.add_config_value('html_theme_options', {}, 'html')
    app.add_config_value('html_title',
                         lambda self: _('%s %s documentation') % (self.project, self.release),
                         'html', [str])
    app.add_config_value('html_short_title', lambda self: self.html_title, 'html')
    app.add_config_value('html_style', None, 'html', [str])
    app.add_config_value('html_logo', None, 'html', [str])
    app.add_config_value('html_favicon', None, 'html', [str])
    app.add_config_value('html_css_files', [], 'html')
    app.add_config_value('html_js_files', [], 'html')
    app.add_config_value('html_static_path', [], 'html')
    app.add_config_value('html_extra_path', [], 'html')
    app.add_config_value('html_last_updated_fmt', None, 'html', [str])
    app.add_config_value('html_sidebars', {}, 'html')
    app.add_config_value('html_additional_pages', {}, 'html')
    app.add_config_value('html_domain_indices', True, 'html', [list])
    app.add_config_value('html_add_permalinks', UNSET, 'html')
    app.add_config_value('html_permalinks', True, 'html')
    app.add_config_value('html_permalinks_icon', '¶', 'html')
    app.add_config_value('html_use_index', True, 'html')
    app.add_config_value('html_split_index', False, 'html')
    app.add_config_value('html_copy_source', True, 'html')
    app.add_config_value('html_show_sourcelink', True, 'html')
    app.add_config_value('html_sourcelink_suffix', '.txt', 'html')
    app.add_config_value('html_use_opensearch', '', 'html')
    app.add_config_value('html_file_suffix', None, 'html', [str])
    app.add_config_value('html_link_suffix', None, 'html', [str])
    app.add_config_value('html_show_copyright', True, 'html')
    app.add_config_value('html_show_sphinx', True, 'html')
    app.add_config_value('html_context', {}, 'html')
    app.add_config_value('html_output_encoding', 'utf-8', 'html')
    app.add_config_value('html_compact_lists', True, 'html')
    app.add_config_value('html_secnumber_suffix', '. ', 'html')
    app.add_config_value('html_search_language', None, 'html', [str])
    app.add_config_value('html_search_options', {}, 'html')
    app.add_config_value('html_search_scorer', '', None)
    app.add_config_value('html_scaled_image_link', True, 'html')
    app.add_config_value('html_baseurl', '', 'html')
    app.add_config_value('html_codeblock_linenos_style', 'inline', 'html',  # RemovedInSphinx60Warning  # NOQA
                         ENUM('table', 'inline'))
    app.add_config_value('html_math_renderer', None, 'env')
    app.add_config_value('html4_writer', False, 'html')

    # events
    app.add_event('html-collect-pages')
    app.add_event('html-page-context')

    # event handlers
    app.connect('config-inited', convert_html_css_files, priority=800)
    app.connect('config-inited', convert_html_js_files, priority=800)
    app.connect('config-inited', migrate_html_add_permalinks, priority=800)
    app.connect('config-inited', validate_html_extra_path, priority=800)
    app.connect('config-inited', validate_html_static_path, priority=800)
    app.connect('config-inited', validate_html_logo, priority=800)
    app.connect('config-inited', validate_html_favicon, priority=800)
    app.connect('builder-inited', validate_math_renderer)
    app.connect('html-page-context', setup_css_tag_helper)
    app.connect('html-page-context', setup_js_tag_helper)
    app.connect('html-page-context', setup_resource_paths)

    # load default math renderer
    app.setup_extension('sphinx.ext.mathjax')

    # load transforms for HTML builder
    app.setup_extension('sphinx.builders.html.transforms')

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
