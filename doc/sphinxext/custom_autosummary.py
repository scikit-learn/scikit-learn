"""
Use custom_autosummary_file_map sphinx configuration option to map filenames
to a custom filename:

custom_autosummary_file_map = {
    "sklearn.cluster.dbscan": "sklearn.cluster.dbscan_lowercase",
    "sklearn.cluster.optics": "sklearn.cluster.optics_lowercase"
}
"""


import os
import logging
from contextlib import suppress

import sphinx
from sphinx import package_dir
from sphinx.jinja2glue import BuiltinTemplateLoader
from jinja2 import FileSystemLoader, TemplateNotFound
from jinja2.sandbox import SandboxedEnvironment

from sphinx.ext.autosummary import Autosummary
from sphinx.ext.autosummary import autosummary_toc
from sphinx.ext.autosummary import autosummary_toc_visit_html
from sphinx.ext.autosummary import autosummary_noop
from sphinx.ext.autosummary import autosummary_table
from sphinx.ext.autosummary import autosummary_table_visit_html
from sphinx.ext.autosummary import autolink_role
from sphinx.ext.autosummary import process_autosummary_toc
from sphinx.ext.autosummary import get_rst_suffix
from sphinx.ext.autosummary import import_by_name
from sphinx.ext.autosummary import get_documenter
from sphinx.ext.autosummary.generate import _underline
from sphinx.ext.autosummary.generate import find_autosummary_in_files
from sphinx.util.inspect import safe_getattr
from sphinx.util.rst import escape as rst_escape
from sphinx.util.osutil import ensuredir

logger = logging.getLogger(__name__)


class CustomAutosummary(Autosummary):
    def get_items(self, names):
        items = super().get_items(names)
        custom_autosummary_file_map = self.config.custom_autosummary_file_map

        output = []
        for name, sig, summary, real_name in items:
            with suppress(KeyError):
                real_name = custom_autosummary_file_map[real_name]
            output.append((name, sig, summary, real_name))

        return output


def generate_autosummary_docs_custom(sources,
                                     output_dir=None,
                                     suffix='.rst',
                                     warn=print,
                                     info=print,
                                     base_path=None,
                                     builder=None,
                                     template_dir=None,
                                     imported_members=False,
                                     app=None):
    showed_sources = list(sorted(sources))
    if len(showed_sources) > 20:
        showed_sources = showed_sources[:10] + ['...'] + showed_sources[-10:]
    info('[autosummary] generating autosummary for: %s' %
         ', '.join(showed_sources))

    if output_dir:
        info('[autosummary] writing to %s' % output_dir)

    if base_path is not None:
        sources = [os.path.join(base_path, filename) for filename in sources]

    # create our own templating environment
    template_dirs = None  # type: List[str]
    template_dirs = [
        os.path.join(package_dir, 'ext', 'autosummary', 'templates')
    ]

    template_loader = None
    if builder is not None:
        # allow the user to override the templates
        template_loader = BuiltinTemplateLoader()
        template_loader.init(builder, dirs=template_dirs)
    else:
        if template_dir:
            template_dirs.insert(0, template_dir)
        template_loader = FileSystemLoader(template_dirs)
    template_env = SandboxedEnvironment(loader=template_loader)
    template_env.filters['underline'] = _underline

    # replace the builtin html filters
    template_env.filters['escape'] = rst_escape
    template_env.filters['e'] = rst_escape

    # read
    items = find_autosummary_in_files(sources)

    # keep track of new files
    new_files = []

    # write
    for name, path, template_name in sorted(set(items), key=str):
        if path is None:
            # The corresponding autosummary:: directive did not have
            # a :toctree: option
            continue

        path = output_dir or os.path.abspath(path)
        ensuredir(path)

        try:
            name, obj, parent, mod_name = import_by_name(name)
        except ImportError as e:
            warn('[autosummary] failed to import %r: %s' % (name, e))
            continue

        try:
            filename = app.config.custom_autosummary_file_map[name] + suffix
        except KeyError:
            filename = name + suffix

        fn = os.path.join(path, filename)

        # skip it if it exists
        if os.path.isfile(fn):
            continue

        new_files.append(fn)

        with open(fn, 'w') as f:
            doc = get_documenter(app, obj, parent)

            if template_name is not None:
                template = template_env.get_template(template_name)
            else:
                try:
                    template = template_env.get_template(
                        'autosummary/%s.rst' % doc.objtype)
                except TemplateNotFound:
                    template = template_env.get_template(
                        'autosummary/base.rst')

            def get_members(obj, typ, include_public=[], imported=True):
                items = []  # type: List[str]
                for name in dir(obj):
                    try:
                        value = safe_getattr(obj, name)
                    except AttributeError:
                        continue
                    documenter = get_documenter(app, value, obj)
                    if documenter.objtype == typ:
                        if imported or getattr(value, '__module__',
                                               None) == obj.__name__:
                            # skip imported members if expected
                            items.append(name)
                public = [
                    x for x in items
                    if x in include_public or not x.startswith('_')
                ]
                return public, items

            ns = {}  # type: Dict[str, Any]

            if doc.objtype == 'module':
                ns['members'] = dir(obj)
                ns['functions'], ns['all_functions'] = \
                    get_members(obj, 'function', imported=imported_members)
                ns['classes'], ns['all_classes'] = \
                    get_members(obj, 'class', imported=imported_members)
                ns['exceptions'], ns['all_exceptions'] = \
                    get_members(obj, 'exception', imported=imported_members)
            elif doc.objtype == 'class':
                ns['members'] = dir(obj)
                ns['inherited_members'] = \
                    set(dir(obj)) - set(obj.__dict__.keys())
                ns['methods'], ns['all_methods'] = \
                    get_members(obj, 'method', ['__init__'])
                ns['attributes'], ns['all_attributes'] = \
                    get_members(obj, 'attribute')

            parts = name.split('.')
            if doc.objtype in ('method', 'attribute'):
                mod_name = '.'.join(parts[:-2])
                cls_name = parts[-2]
                obj_name = '.'.join(parts[-2:])
                ns['class'] = cls_name
            else:
                mod_name, obj_name = '.'.join(parts[:-1]), parts[-1]

            ns['fullname'] = name
            ns['module'] = mod_name
            ns['objname'] = obj_name
            ns['name'] = parts[-1]

            ns['objtype'] = doc.objtype
            ns['underline'] = len(name) * '='

            rendered = template.render(**ns)
            f.write(rendered)

    # descend recursively to new files
    if new_files:
        generate_autosummary_docs_custom(
            new_files,
            output_dir=output_dir,
            suffix=suffix,
            warn=warn,
            info=info,
            base_path=base_path,
            builder=builder,
            template_dir=template_dir,
            app=app)


def process_generate_options_custom(app):
    # type: (Sphinx) -> None
    genfiles = app.config.autosummary_generate

    if genfiles and not hasattr(genfiles, '__len__'):
        env = app.builder.env
        genfiles = [
            env.doc2path(x, base=None) for x in env.found_docs
            if os.path.isfile(env.doc2path(x))
        ]

    if not genfiles:
        return

    ext = list(app.config.source_suffix)
    genfiles = [
        genfile + (not genfile.endswith(tuple(ext)) and ext[0] or '')
        for genfile in genfiles
    ]

    suffix = get_rst_suffix(app)
    if suffix is None:
        logger.warning(
            ('autosummary generats .rst files internally. '
             'But your source_suffix does not contain .rst. Skipped.'))
        return

    generate_autosummary_docs_custom(
        genfiles,
        builder=app.builder,
        warn=logger.warning,
        info=logger.info,
        suffix=suffix,
        base_path=app.srcdir,
        app=app)


def setup(app):
    app.setup_extension('sphinx.ext.autodoc')
    app.add_config_value("custom_autosummary_file_map", {}, None)

    app.add_node(
        autosummary_toc,
        html=(autosummary_toc_visit_html, autosummary_noop),
        latex=(autosummary_noop, autosummary_noop),
        text=(autosummary_noop, autosummary_noop),
        man=(autosummary_noop, autosummary_noop),
        texinfo=(autosummary_noop, autosummary_noop))
    app.add_node(
        autosummary_table,
        html=(autosummary_table_visit_html, autosummary_noop),
        latex=(autosummary_noop, autosummary_noop),
        text=(autosummary_noop, autosummary_noop),
        man=(autosummary_noop, autosummary_noop),
        texinfo=(autosummary_noop, autosummary_noop))
    app.add_directive("autosummary", CustomAutosummary)
    app.add_role('autolink', autolink_role)
    app.connect('doctree-read', process_autosummary_toc)
    app.connect('builder-inited', process_generate_options_custom)
    app.add_config_value('autosummary_generate', [], True, [bool])
    return {'version': sphinx.__display_version__, 'parallel_read_safe': True}
