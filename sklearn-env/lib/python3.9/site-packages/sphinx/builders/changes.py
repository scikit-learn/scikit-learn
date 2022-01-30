"""
    sphinx.builders.changes
    ~~~~~~~~~~~~~~~~~~~~~~~

    Changelog builder.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import html
from os import path
from typing import Any, Dict, List, Tuple, cast

from sphinx import package_dir
from sphinx.application import Sphinx
from sphinx.builders import Builder
from sphinx.domains.changeset import ChangeSetDomain
from sphinx.locale import _, __
from sphinx.theming import HTMLThemeFactory
from sphinx.util import logging
from sphinx.util.console import bold  # type: ignore
from sphinx.util.fileutil import copy_asset_file
from sphinx.util.osutil import ensuredir, os_path

logger = logging.getLogger(__name__)


class ChangesBuilder(Builder):
    """
    Write a summary with all versionadded/changed directives.
    """
    name = 'changes'
    epilog = __('The overview file is in %(outdir)s.')

    def init(self) -> None:
        self.create_template_bridge()
        theme_factory = HTMLThemeFactory(self.app)
        self.theme = theme_factory.create('default')
        self.templates.init(self, self.theme)

    def get_outdated_docs(self) -> str:
        return self.outdir

    typemap = {
        'versionadded': 'added',
        'versionchanged': 'changed',
        'deprecated': 'deprecated',
    }

    def write(self, *ignored: Any) -> None:
        version = self.config.version
        domain = cast(ChangeSetDomain, self.env.get_domain('changeset'))
        libchanges: Dict[str, List[Tuple[str, str, int]]] = {}
        apichanges: List[Tuple[str, str, int]] = []
        otherchanges: Dict[Tuple[str, str], List[Tuple[str, str, int]]] = {}

        changesets = domain.get_changesets_for(version)
        if not changesets:
            logger.info(bold(__('no changes in version %s.') % version))
            return
        logger.info(bold(__('writing summary file...')))
        for changeset in changesets:
            if isinstance(changeset.descname, tuple):
                descname = changeset.descname[0]
            else:
                descname = changeset.descname
            ttext = self.typemap[changeset.type]
            context = changeset.content.replace('\n', ' ')
            if descname and changeset.docname.startswith('c-api'):
                if context:
                    entry = '<b>%s</b>: <i>%s:</i> %s' % (descname, ttext,
                                                          context)
                else:
                    entry = '<b>%s</b>: <i>%s</i>.' % (descname, ttext)
                apichanges.append((entry, changeset.docname, changeset.lineno))
            elif descname or changeset.module:
                module = changeset.module or _('Builtins')
                if not descname:
                    descname = _('Module level')
                if context:
                    entry = '<b>%s</b>: <i>%s:</i> %s' % (descname, ttext,
                                                          context)
                else:
                    entry = '<b>%s</b>: <i>%s</i>.' % (descname, ttext)
                libchanges.setdefault(module, []).append((entry, changeset.docname,
                                                          changeset.lineno))
            else:
                if not context:
                    continue
                entry = '<i>%s:</i> %s' % (ttext.capitalize(), context)
                title = self.env.titles[changeset.docname].astext()
                otherchanges.setdefault((changeset.docname, title), []).append(
                    (entry, changeset.docname, changeset.lineno))

        ctx = {
            'project': self.config.project,
            'version': version,
            'docstitle': self.config.html_title,
            'shorttitle': self.config.html_short_title,
            'libchanges': sorted(libchanges.items()),
            'apichanges': sorted(apichanges),
            'otherchanges': sorted(otherchanges.items()),
            'show_copyright': self.config.html_show_copyright,
            'show_sphinx': self.config.html_show_sphinx,
        }
        with open(path.join(self.outdir, 'index.html'), 'w', encoding='utf8') as f:
            f.write(self.templates.render('changes/frameset.html', ctx))
        with open(path.join(self.outdir, 'changes.html'), 'w', encoding='utf8') as f:
            f.write(self.templates.render('changes/versionchanges.html', ctx))

        hltext = ['.. versionadded:: %s' % version,
                  '.. versionchanged:: %s' % version,
                  '.. deprecated:: %s' % version]

        def hl(no: int, line: str) -> str:
            line = '<a name="L%s"> </a>' % no + html.escape(line)
            for x in hltext:
                if x in line:
                    line = '<span class="hl">%s</span>' % line
                    break
            return line

        logger.info(bold(__('copying source files...')))
        for docname in self.env.all_docs:
            with open(self.env.doc2path(docname),
                      encoding=self.env.config.source_encoding) as f:
                try:
                    lines = f.readlines()
                except UnicodeDecodeError:
                    logger.warning(__('could not read %r for changelog creation'), docname)
                    continue
            targetfn = path.join(self.outdir, 'rst', os_path(docname)) + '.html'
            ensuredir(path.dirname(targetfn))
            with open(targetfn, 'w', encoding='utf-8') as f:
                text = ''.join(hl(i + 1, line) for (i, line) in enumerate(lines))
                ctx = {
                    'filename': self.env.doc2path(docname, None),
                    'text': text
                }
                f.write(self.templates.render('changes/rstsource.html', ctx))
        themectx = {'theme_' + key: val for (key, val) in
                    self.theme.get_options({}).items()}
        copy_asset_file(path.join(package_dir, 'themes', 'default', 'static', 'default.css_t'),
                        self.outdir, context=themectx, renderer=self.templates)
        copy_asset_file(path.join(package_dir, 'themes', 'basic', 'static', 'basic.css'),
                        self.outdir)

    def hl(self, text: str, version: str) -> str:
        text = html.escape(text)
        for directive in ('versionchanged', 'versionadded', 'deprecated'):
            text = text.replace('.. %s:: %s' % (directive, version),
                                '<b>.. %s:: %s</b>' % (directive, version))
        return text

    def finish(self) -> None:
        pass


def setup(app: Sphinx) -> Dict[str, Any]:
    app.add_builder(ChangesBuilder)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
