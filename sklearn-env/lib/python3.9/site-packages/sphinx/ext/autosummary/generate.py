"""
    sphinx.ext.autosummary.generate
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Usable as a library or script to generate automatic RST source files for
    items referred to in autosummary:: directives.

    Each generated RST file contains a single auto*:: directive which
    extracts the docstring of the referred item.

    Example Makefile rule::

       generate:
               sphinx-autogen -o source/generated source/*.rst

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import argparse
import inspect
import locale
import os
import pkgutil
import pydoc
import re
import sys
import warnings
from gettext import NullTranslations
from os import path
from typing import Any, Dict, List, NamedTuple, Sequence, Set, Tuple, Type, Union

from jinja2 import TemplateNotFound
from jinja2.sandbox import SandboxedEnvironment

import sphinx.locale
from sphinx import __display_version__, package_dir
from sphinx.application import Sphinx
from sphinx.builders import Builder
from sphinx.config import Config
from sphinx.deprecation import RemovedInSphinx50Warning
from sphinx.ext.autodoc import Documenter
from sphinx.ext.autodoc.importer import import_module
from sphinx.ext.autosummary import (ImportExceptionGroup, get_documenter, import_by_name,
                                    import_ivar_by_name)
from sphinx.locale import __
from sphinx.pycode import ModuleAnalyzer, PycodeError
from sphinx.registry import SphinxComponentRegistry
from sphinx.util import logging, rst, split_full_qualified_name
from sphinx.util.inspect import getall, safe_getattr
from sphinx.util.osutil import ensuredir
from sphinx.util.template import SphinxTemplateLoader

logger = logging.getLogger(__name__)


class DummyApplication:
    """Dummy Application class for sphinx-autogen command."""

    def __init__(self, translator: NullTranslations) -> None:
        self.config = Config()
        self.registry = SphinxComponentRegistry()
        self.messagelog: List[str] = []
        self.srcdir = "/"
        self.translator = translator
        self.verbosity = 0
        self._warncount = 0
        self.warningiserror = False

        self.config.add('autosummary_context', {}, True, None)
        self.config.add('autosummary_filename_map', {}, True, None)
        self.config.add('autosummary_ignore_module_all', True, 'env', bool)
        self.config.init_values()

    def emit_firstresult(self, *args: Any) -> None:
        pass


class AutosummaryEntry(NamedTuple):
    name: str
    path: str
    template: str
    recursive: bool


def setup_documenters(app: Any) -> None:
    from sphinx.ext.autodoc import (AttributeDocumenter, ClassDocumenter, DataDocumenter,
                                    DecoratorDocumenter, ExceptionDocumenter,
                                    FunctionDocumenter, MethodDocumenter, ModuleDocumenter,
                                    NewTypeAttributeDocumenter, NewTypeDataDocumenter,
                                    PropertyDocumenter)
    documenters: List[Type[Documenter]] = [
        ModuleDocumenter, ClassDocumenter, ExceptionDocumenter, DataDocumenter,
        FunctionDocumenter, MethodDocumenter, NewTypeAttributeDocumenter,
        NewTypeDataDocumenter, AttributeDocumenter, DecoratorDocumenter, PropertyDocumenter,
    ]
    for documenter in documenters:
        app.registry.add_documenter(documenter.objtype, documenter)


def _simple_info(msg: str) -> None:
    warnings.warn('_simple_info() is deprecated.',
                  RemovedInSphinx50Warning, stacklevel=2)
    print(msg)


def _simple_warn(msg: str) -> None:
    warnings.warn('_simple_warn() is deprecated.',
                  RemovedInSphinx50Warning, stacklevel=2)
    print('WARNING: ' + msg, file=sys.stderr)


def _underline(title: str, line: str = '=') -> str:
    if '\n' in title:
        raise ValueError('Can only underline single lines')
    return title + '\n' + line * len(title)


class AutosummaryRenderer:
    """A helper class for rendering."""

    def __init__(self, app: Union[Builder, Sphinx], template_dir: str = None) -> None:
        if isinstance(app, Builder):
            warnings.warn('The first argument for AutosummaryRenderer has been '
                          'changed to Sphinx object',
                          RemovedInSphinx50Warning, stacklevel=2)
        if template_dir:
            warnings.warn('template_dir argument for AutosummaryRenderer is deprecated.',
                          RemovedInSphinx50Warning, stacklevel=2)

        system_templates_path = [os.path.join(package_dir, 'ext', 'autosummary', 'templates')]
        loader = SphinxTemplateLoader(app.srcdir, app.config.templates_path,
                                      system_templates_path)

        self.env = SandboxedEnvironment(loader=loader)
        self.env.filters['escape'] = rst.escape
        self.env.filters['e'] = rst.escape
        self.env.filters['underline'] = _underline

        if isinstance(app, (Sphinx, DummyApplication)):
            if app.translator:
                self.env.add_extension("jinja2.ext.i18n")
                self.env.install_gettext_translations(app.translator)
        elif isinstance(app, Builder):
            if app.app.translator:
                self.env.add_extension("jinja2.ext.i18n")
                self.env.install_gettext_translations(app.app.translator)

    def exists(self, template_name: str) -> bool:
        """Check if template file exists."""
        warnings.warn('AutosummaryRenderer.exists() is deprecated.',
                      RemovedInSphinx50Warning, stacklevel=2)
        try:
            self.env.get_template(template_name)
            return True
        except TemplateNotFound:
            return False

    def render(self, template_name: str, context: Dict) -> str:
        """Render a template file."""
        try:
            template = self.env.get_template(template_name)
        except TemplateNotFound:
            try:
                # objtype is given as template_name
                template = self.env.get_template('autosummary/%s.rst' % template_name)
            except TemplateNotFound:
                # fallback to base.rst
                template = self.env.get_template('autosummary/base.rst')

        return template.render(context)


# -- Generating output ---------------------------------------------------------


class ModuleScanner:
    def __init__(self, app: Any, obj: Any) -> None:
        self.app = app
        self.object = obj

    def get_object_type(self, name: str, value: Any) -> str:
        return get_documenter(self.app, value, self.object).objtype

    def is_skipped(self, name: str, value: Any, objtype: str) -> bool:
        try:
            return self.app.emit_firstresult('autodoc-skip-member', objtype,
                                             name, value, False, {})
        except Exception as exc:
            logger.warning(__('autosummary: failed to determine %r to be documented, '
                              'the following exception was raised:\n%s'),
                           name, exc, type='autosummary')
            return False

    def scan(self, imported_members: bool) -> List[str]:
        members = []
        for name in members_of(self.object, self.app.config):
            try:
                value = safe_getattr(self.object, name)
            except AttributeError:
                value = None

            objtype = self.get_object_type(name, value)
            if self.is_skipped(name, value, objtype):
                continue

            try:
                if inspect.ismodule(value):
                    imported = True
                elif safe_getattr(value, '__module__') != self.object.__name__:
                    imported = True
                else:
                    imported = False
            except AttributeError:
                imported = False

            respect_module_all = not self.app.config.autosummary_ignore_module_all
            if imported_members:
                # list all members up
                members.append(name)
            elif imported is False:
                # list not-imported members
                members.append(name)
            elif '__all__' in dir(self.object) and respect_module_all:
                # list members that have __all__ set
                members.append(name)

        return members


def members_of(obj: Any, conf: Config) -> Sequence[str]:
    """Get the members of ``obj``, possibly ignoring the ``__all__`` module attribute

    Follows the ``conf.autosummary_ignore_module_all`` setting."""

    if conf.autosummary_ignore_module_all:
        return dir(obj)
    else:
        return getall(obj) or dir(obj)


def generate_autosummary_content(name: str, obj: Any, parent: Any,
                                 template: AutosummaryRenderer, template_name: str,
                                 imported_members: bool, app: Any,
                                 recursive: bool, context: Dict,
                                 modname: str = None, qualname: str = None) -> str:
    doc = get_documenter(app, obj, parent)

    def skip_member(obj: Any, name: str, objtype: str) -> bool:
        try:
            return app.emit_firstresult('autodoc-skip-member', objtype, name,
                                        obj, False, {})
        except Exception as exc:
            logger.warning(__('autosummary: failed to determine %r to be documented, '
                              'the following exception was raised:\n%s'),
                           name, exc, type='autosummary')
            return False

    def get_class_members(obj: Any) -> Dict[str, Any]:
        members = sphinx.ext.autodoc.get_class_members(obj, [qualname], safe_getattr)
        return {name: member.object for name, member in members.items()}

    def get_module_members(obj: Any) -> Dict[str, Any]:
        members = {}
        for name in members_of(obj, app.config):
            try:
                members[name] = safe_getattr(obj, name)
            except AttributeError:
                continue
        return members

    def get_all_members(obj: Any) -> Dict[str, Any]:
        if doc.objtype == "module":
            return get_module_members(obj)
        elif doc.objtype == "class":
            return get_class_members(obj)
        return {}

    def get_members(obj: Any, types: Set[str], include_public: List[str] = [],
                    imported: bool = True) -> Tuple[List[str], List[str]]:
        items: List[str] = []
        public: List[str] = []

        all_members = get_all_members(obj)
        for name, value in all_members.items():
            documenter = get_documenter(app, value, obj)
            if documenter.objtype in types:
                # skip imported members if expected
                if imported or getattr(value, '__module__', None) == obj.__name__:
                    skipped = skip_member(value, name, documenter.objtype)
                    if skipped is True:
                        pass
                    elif skipped is False:
                        # show the member forcedly
                        items.append(name)
                        public.append(name)
                    else:
                        items.append(name)
                        if name in include_public or not name.startswith('_'):
                            # considers member as public
                            public.append(name)
        return public, items

    def get_module_attrs(members: Any) -> Tuple[List[str], List[str]]:
        """Find module attributes with docstrings."""
        attrs, public = [], []
        try:
            analyzer = ModuleAnalyzer.for_module(name)
            attr_docs = analyzer.find_attr_docs()
            for namespace, attr_name in attr_docs:
                if namespace == '' and attr_name in members:
                    attrs.append(attr_name)
                    if not attr_name.startswith('_'):
                        public.append(attr_name)
        except PycodeError:
            pass    # give up if ModuleAnalyzer fails to parse code
        return public, attrs

    def get_modules(obj: Any) -> Tuple[List[str], List[str]]:
        items: List[str] = []
        for _, modname, _ispkg in pkgutil.iter_modules(obj.__path__):
            fullname = name + '.' + modname
            try:
                module = import_module(fullname)
                if module and hasattr(module, '__sphinx_mock__'):
                    continue
            except ImportError:
                pass

            items.append(fullname)
        public = [x for x in items if not x.split('.')[-1].startswith('_')]
        return public, items

    ns: Dict[str, Any] = {}
    ns.update(context)

    if doc.objtype == 'module':
        scanner = ModuleScanner(app, obj)
        ns['members'] = scanner.scan(imported_members)
        ns['functions'], ns['all_functions'] = \
            get_members(obj, {'function'}, imported=imported_members)
        ns['classes'], ns['all_classes'] = \
            get_members(obj, {'class'}, imported=imported_members)
        ns['exceptions'], ns['all_exceptions'] = \
            get_members(obj, {'exception'}, imported=imported_members)
        ns['attributes'], ns['all_attributes'] = \
            get_module_attrs(ns['members'])
        ispackage = hasattr(obj, '__path__')
        if ispackage and recursive:
            ns['modules'], ns['all_modules'] = get_modules(obj)
    elif doc.objtype == 'class':
        ns['members'] = dir(obj)
        ns['inherited_members'] = \
            set(dir(obj)) - set(obj.__dict__.keys())
        ns['methods'], ns['all_methods'] = \
            get_members(obj, {'method'}, ['__init__'])
        ns['attributes'], ns['all_attributes'] = \
            get_members(obj, {'attribute', 'property'})

    if modname is None or qualname is None:
        modname, qualname = split_full_qualified_name(name)

    if doc.objtype in ('method', 'attribute', 'property'):
        ns['class'] = qualname.rsplit(".", 1)[0]

    if doc.objtype in ('class',):
        shortname = qualname
    else:
        shortname = qualname.rsplit(".", 1)[-1]

    ns['fullname'] = name
    ns['module'] = modname
    ns['objname'] = qualname
    ns['name'] = shortname

    ns['objtype'] = doc.objtype
    ns['underline'] = len(name) * '='

    if template_name:
        return template.render(template_name, ns)
    else:
        return template.render(doc.objtype, ns)


def generate_autosummary_docs(sources: List[str], output_dir: str = None,
                              suffix: str = '.rst', base_path: str = None,
                              builder: Builder = None, template_dir: str = None,
                              imported_members: bool = False, app: Any = None,
                              overwrite: bool = True, encoding: str = 'utf-8') -> None:
    if builder:
        warnings.warn('builder argument for generate_autosummary_docs() is deprecated.',
                      RemovedInSphinx50Warning, stacklevel=2)

    if template_dir:
        warnings.warn('template_dir argument for generate_autosummary_docs() is deprecated.',
                      RemovedInSphinx50Warning, stacklevel=2)

    showed_sources = list(sorted(sources))
    if len(showed_sources) > 20:
        showed_sources = showed_sources[:10] + ['...'] + showed_sources[-10:]
    logger.info(__('[autosummary] generating autosummary for: %s') %
                ', '.join(showed_sources))

    if output_dir:
        logger.info(__('[autosummary] writing to %s') % output_dir)

    if base_path is not None:
        sources = [os.path.join(base_path, filename) for filename in sources]

    template = AutosummaryRenderer(app)

    # read
    items = find_autosummary_in_files(sources)

    # keep track of new files
    new_files = []

    if app:
        filename_map = app.config.autosummary_filename_map
    else:
        filename_map = {}

    # write
    for entry in sorted(set(items), key=str):
        if entry.path is None:
            # The corresponding autosummary:: directive did not have
            # a :toctree: option
            continue

        path = output_dir or os.path.abspath(entry.path)
        ensuredir(path)

        try:
            name, obj, parent, modname = import_by_name(entry.name, grouped_exception=True)
            qualname = name.replace(modname + ".", "")
        except ImportExceptionGroup as exc:
            try:
                # try to import as an instance attribute
                name, obj, parent, modname = import_ivar_by_name(entry.name)
                qualname = name.replace(modname + ".", "")
            except ImportError as exc2:
                if exc2.__cause__:
                    exceptions: List[BaseException] = exc.exceptions + [exc2.__cause__]
                else:
                    exceptions = exc.exceptions + [exc2]

                errors = list(set("* %s: %s" % (type(e).__name__, e) for e in exceptions))
                logger.warning(__('[autosummary] failed to import %s.\nPossible hints:\n%s'),
                               entry.name, '\n'.join(errors))
                continue

        context: Dict[str, Any] = {}
        if app:
            context.update(app.config.autosummary_context)

        content = generate_autosummary_content(name, obj, parent, template, entry.template,
                                               imported_members, app, entry.recursive, context,
                                               modname, qualname)

        filename = os.path.join(path, filename_map.get(name, name) + suffix)
        if os.path.isfile(filename):
            with open(filename, encoding=encoding) as f:
                old_content = f.read()

            if content == old_content:
                continue
            elif overwrite:  # content has changed
                with open(filename, 'w', encoding=encoding) as f:
                    f.write(content)
                new_files.append(filename)
        else:
            with open(filename, 'w', encoding=encoding) as f:
                f.write(content)
            new_files.append(filename)

    # descend recursively to new files
    if new_files:
        generate_autosummary_docs(new_files, output_dir=output_dir,
                                  suffix=suffix, base_path=base_path,
                                  builder=builder, template_dir=template_dir,
                                  imported_members=imported_members, app=app,
                                  overwrite=overwrite)


# -- Finding documented entries in files ---------------------------------------

def find_autosummary_in_files(filenames: List[str]) -> List[AutosummaryEntry]:
    """Find out what items are documented in source/*.rst.

    See `find_autosummary_in_lines`.
    """
    documented: List[AutosummaryEntry] = []
    for filename in filenames:
        with open(filename, encoding='utf-8', errors='ignore') as f:
            lines = f.read().splitlines()
            documented.extend(find_autosummary_in_lines(lines, filename=filename))
    return documented


def find_autosummary_in_docstring(name: str, module: str = None, filename: str = None
                                  ) -> List[AutosummaryEntry]:
    """Find out what items are documented in the given object's docstring.

    See `find_autosummary_in_lines`.
    """
    if module:
        warnings.warn('module argument for find_autosummary_in_docstring() is deprecated.',
                      RemovedInSphinx50Warning, stacklevel=2)

    try:
        real_name, obj, parent, modname = import_by_name(name, grouped_exception=True)
        lines = pydoc.getdoc(obj).splitlines()
        return find_autosummary_in_lines(lines, module=name, filename=filename)
    except AttributeError:
        pass
    except ImportExceptionGroup as exc:
        errors = list(set("* %s: %s" % (type(e).__name__, e) for e in exc.exceptions))
        print('Failed to import %s.\nPossible hints:\n%s' % (name, '\n'.join(errors)))
    except SystemExit:
        print("Failed to import '%s'; the module executes module level "
              "statement and it might call sys.exit()." % name)
    return []


def find_autosummary_in_lines(lines: List[str], module: str = None, filename: str = None
                              ) -> List[AutosummaryEntry]:
    """Find out what items appear in autosummary:: directives in the
    given lines.

    Returns a list of (name, toctree, template) where *name* is a name
    of an object and *toctree* the :toctree: path of the corresponding
    autosummary directive (relative to the root of the file name), and
    *template* the value of the :template: option. *toctree* and
    *template* ``None`` if the directive does not have the
    corresponding options set.
    """
    autosummary_re = re.compile(r'^(\s*)\.\.\s+autosummary::\s*')
    automodule_re = re.compile(
        r'^\s*\.\.\s+automodule::\s*([A-Za-z0-9_.]+)\s*$')
    module_re = re.compile(
        r'^\s*\.\.\s+(current)?module::\s*([a-zA-Z0-9_.]+)\s*$')
    autosummary_item_re = re.compile(r'^\s+(~?[_a-zA-Z][a-zA-Z0-9_.]*)\s*.*?')
    recursive_arg_re = re.compile(r'^\s+:recursive:\s*$')
    toctree_arg_re = re.compile(r'^\s+:toctree:\s*(.*?)\s*$')
    template_arg_re = re.compile(r'^\s+:template:\s*(.*?)\s*$')

    documented: List[AutosummaryEntry] = []

    recursive = False
    toctree: str = None
    template = None
    current_module = module
    in_autosummary = False
    base_indent = ""

    for line in lines:
        if in_autosummary:
            m = recursive_arg_re.match(line)
            if m:
                recursive = True
                continue

            m = toctree_arg_re.match(line)
            if m:
                toctree = m.group(1)
                if filename:
                    toctree = os.path.join(os.path.dirname(filename),
                                           toctree)
                continue

            m = template_arg_re.match(line)
            if m:
                template = m.group(1).strip()
                continue

            if line.strip().startswith(':'):
                continue  # skip options

            m = autosummary_item_re.match(line)
            if m:
                name = m.group(1).strip()
                if name.startswith('~'):
                    name = name[1:]
                if current_module and \
                   not name.startswith(current_module + '.'):
                    name = "%s.%s" % (current_module, name)
                documented.append(AutosummaryEntry(name, toctree, template, recursive))
                continue

            if not line.strip() or line.startswith(base_indent + " "):
                continue

            in_autosummary = False

        m = autosummary_re.match(line)
        if m:
            in_autosummary = True
            base_indent = m.group(1)
            recursive = False
            toctree = None
            template = None
            continue

        m = automodule_re.search(line)
        if m:
            current_module = m.group(1).strip()
            # recurse into the automodule docstring
            documented.extend(find_autosummary_in_docstring(
                current_module, filename=filename))
            continue

        m = module_re.match(line)
        if m:
            current_module = m.group(2)
            continue

    return documented


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage='%(prog)s [OPTIONS] <SOURCE_FILE>...',
        epilog=__('For more information, visit <https://www.sphinx-doc.org/>.'),
        description=__("""
Generate ReStructuredText using autosummary directives.

sphinx-autogen is a frontend to sphinx.ext.autosummary.generate. It generates
the reStructuredText files from the autosummary directives contained in the
given input files.

The format of the autosummary directive is documented in the
``sphinx.ext.autosummary`` Python module and can be read using::

  pydoc sphinx.ext.autosummary
"""))

    parser.add_argument('--version', action='version', dest='show_version',
                        version='%%(prog)s %s' % __display_version__)

    parser.add_argument('source_file', nargs='+',
                        help=__('source files to generate rST files for'))

    parser.add_argument('-o', '--output-dir', action='store',
                        dest='output_dir',
                        help=__('directory to place all output in'))
    parser.add_argument('-s', '--suffix', action='store', dest='suffix',
                        default='rst',
                        help=__('default suffix for files (default: '
                                '%(default)s)'))
    parser.add_argument('-t', '--templates', action='store', dest='templates',
                        default=None,
                        help=__('custom template directory (default: '
                                '%(default)s)'))
    parser.add_argument('-i', '--imported-members', action='store_true',
                        dest='imported_members', default=False,
                        help=__('document imported members (default: '
                                '%(default)s)'))
    parser.add_argument('-a', '--respect-module-all', action='store_true',
                        dest='respect_module_all', default=False,
                        help=__('document exactly the members in module __all__ attribute. '
                                '(default: %(default)s)'))

    return parser


def main(argv: List[str] = sys.argv[1:]) -> None:
    sphinx.locale.setlocale(locale.LC_ALL, '')
    sphinx.locale.init_console(os.path.join(package_dir, 'locale'), 'sphinx')
    translator, _ = sphinx.locale.init([], None)

    app = DummyApplication(translator)
    logging.setup(app, sys.stdout, sys.stderr)  # type: ignore
    setup_documenters(app)
    args = get_parser().parse_args(argv)

    if args.templates:
        app.config.templates_path.append(path.abspath(args.templates))
    app.config.autosummary_ignore_module_all = not args.respect_module_all  # type: ignore

    generate_autosummary_docs(args.source_file, args.output_dir,
                              '.' + args.suffix,
                              imported_members=args.imported_members,
                              app=app)


if __name__ == '__main__':
    main()
