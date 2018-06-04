# -*- coding: utf-8 -*-
"""
    sphinx.registry
    ~~~~~~~~~~~~~~~

    Sphinx component registry.

    :copyright: Copyright 2007-2016 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""
from __future__ import print_function

import traceback

from pkg_resources import iter_entry_points
from six import iteritems, itervalues, string_types

from sphinx.domains import ObjType
from sphinx.domains.std import GenericObject, Target
from sphinx.errors import ExtensionError, SphinxError, VersionRequirementError
from sphinx.extension import Extension
from sphinx.locale import __
from sphinx.parsers import Parser as SphinxParser
from sphinx.roles import XRefRole
from sphinx.util import import_object
from sphinx.util import logging
from sphinx.util.console import bold  # type: ignore
from sphinx.util.docutils import directive_helper

if False:
    # For type annotation
    from typing import Any, Callable, Dict, Iterator, List, Type, Union  # NOQA
    from docutils import nodes  # NOQA
    from docutils.io import Input  # NOQA
    from docutils.parsers import Parser  # NOQA
    from docutils.transforms import Transform  # NOQA
    from sphinx.application import Sphinx  # NOQA
    from sphinx.builders import Builder  # NOQA
    from sphinx.domains import Domain, Index  # NOQA
    from sphinx.environment import BuildEnvironment  # NOQA
    from sphinx.ext.autodoc import Documenter  # NOQA
    from sphinx.util.typing import RoleFunction  # NOQA

logger = logging.getLogger(__name__)

# list of deprecated extensions. Keys are extension name.
# Values are Sphinx version that merge the extension.
EXTENSION_BLACKLIST = {
    "sphinxjp.themecore": "1.2"
}  # type: Dict[unicode, unicode]


class SphinxComponentRegistry(object):
    def __init__(self):
        self.autodoc_attrgettrs = {}    # type: Dict[Type, Callable[[Any, unicode, Any], Any]]
        self.builders = {}              # type: Dict[unicode, Type[Builder]]
        self.documenters = {}           # type: Dict[unicode, Type[Documenter]]
        self.domains = {}               # type: Dict[unicode, Type[Domain]]
        self.domain_directives = {}     # type: Dict[unicode, Dict[unicode, Any]]
        self.domain_indices = {}        # type: Dict[unicode, List[Type[Index]]]
        self.domain_object_types = {}   # type: Dict[unicode, Dict[unicode, ObjType]]
        self.domain_roles = {}          # type: Dict[unicode, Dict[unicode, Union[RoleFunction, XRefRole]]]  # NOQA
        self.post_transforms = []       # type: List[Type[Transform]]
        self.source_parsers = {}        # type: Dict[unicode, Parser]
        self.source_inputs = {}         # type: Dict[unicode, Input]
        self.translators = {}           # type: Dict[unicode, nodes.NodeVisitor]
        self.transforms = []            # type: List[Type[Transform]]

    def add_builder(self, builder):
        # type: (Type[Builder]) -> None
        logger.debug('[app] adding builder: %r', builder)
        if not hasattr(builder, 'name'):
            raise ExtensionError(__('Builder class %s has no "name" attribute') % builder)
        if builder.name in self.builders:
            raise ExtensionError(__('Builder %r already exists (in module %s)') %
                                 (builder.name, self.builders[builder.name].__module__))
        self.builders[builder.name] = builder

    def preload_builder(self, app, name):
        # type: (Sphinx, unicode) -> None
        if name is None:
            return

        if name not in self.builders:
            entry_points = iter_entry_points('sphinx.builders', name)
            try:
                entry_point = next(entry_points)
            except StopIteration:
                raise SphinxError(__('Builder name %s not registered or available'
                                     ' through entry point') % name)

            self.load_extension(app, entry_point.module_name)

    def create_builder(self, app, name):
        # type: (Sphinx, unicode) -> Builder
        if name not in self.builders:
            raise SphinxError(__('Builder name %s not registered') % name)

        return self.builders[name](app)

    def add_domain(self, domain):
        # type: (Type[Domain]) -> None
        logger.debug('[app] adding domain: %r', domain)
        if domain.name in self.domains:
            raise ExtensionError(__('domain %s already registered') % domain.name)
        self.domains[domain.name] = domain

    def has_domain(self, domain):
        # type: (unicode) -> bool
        return domain in self.domains

    def create_domains(self, env):
        # type: (BuildEnvironment) -> Iterator[Domain]
        for DomainClass in itervalues(self.domains):
            domain = DomainClass(env)

            # transplant components added by extensions
            domain.directives.update(self.domain_directives.get(domain.name, {}))
            domain.roles.update(self.domain_roles.get(domain.name, {}))
            domain.indices.extend(self.domain_indices.get(domain.name, []))
            for name, objtype in iteritems(self.domain_object_types.get(domain.name, {})):
                domain.add_object_type(name, objtype)

            yield domain

    def override_domain(self, domain):
        # type: (Type[Domain]) -> None
        logger.debug('[app] overriding domain: %r', domain)
        if domain.name not in self.domains:
            raise ExtensionError(__('domain %s not yet registered') % domain.name)
        if not issubclass(domain, self.domains[domain.name]):
            raise ExtensionError(__('new domain not a subclass of registered %s '
                                    'domain') % domain.name)
        self.domains[domain.name] = domain

    def add_directive_to_domain(self, domain, name, obj,
                                has_content=None, argument_spec=None, **option_spec):
        # type: (unicode, unicode, Any, bool, Any, Any) -> None
        logger.debug('[app] adding directive to domain: %r',
                     (domain, name, obj, has_content, argument_spec, option_spec))
        if domain not in self.domains:
            raise ExtensionError(__('domain %s not yet registered') % domain)
        directives = self.domain_directives.setdefault(domain, {})
        directives[name] = directive_helper(obj, has_content, argument_spec, **option_spec)

    def add_role_to_domain(self, domain, name, role):
        # type: (unicode, unicode, Union[RoleFunction, XRefRole]) -> None
        logger.debug('[app] adding role to domain: %r', (domain, name, role))
        if domain not in self.domains:
            raise ExtensionError(__('domain %s not yet registered') % domain)
        roles = self.domain_roles.setdefault(domain, {})
        roles[name] = role

    def add_index_to_domain(self, domain, index):
        # type: (unicode, Type[Index]) -> None
        logger.debug('[app] adding index to domain: %r', (domain, index))
        if domain not in self.domains:
            raise ExtensionError(__('domain %s not yet registered') % domain)
        indices = self.domain_indices.setdefault(domain, [])
        indices.append(index)

    def add_object_type(self, directivename, rolename, indextemplate='',
                        parse_node=None, ref_nodeclass=None, objname='',
                        doc_field_types=[]):
        # type: (unicode, unicode, unicode, Callable, nodes.Node, unicode, List) -> None
        logger.debug('[app] adding object type: %r',
                     (directivename, rolename, indextemplate, parse_node,
                      ref_nodeclass, objname, doc_field_types))

        # create a subclass of GenericObject as the new directive
        directive = type(directivename,  # type: ignore
                         (GenericObject, object),
                         {'indextemplate': indextemplate,
                          'parse_node': staticmethod(parse_node),
                          'doc_field_types': doc_field_types})

        self.add_directive_to_domain('std', directivename, directive)
        self.add_role_to_domain('std', rolename, XRefRole(innernodeclass=ref_nodeclass))

        object_types = self.domain_object_types.setdefault('std', {})
        object_types[directivename] = ObjType(objname or directivename, rolename)

    def add_crossref_type(self, directivename, rolename, indextemplate='',
                          ref_nodeclass=None, objname=''):
        # type: (unicode, unicode, unicode, nodes.Node, unicode) -> None
        logger.debug('[app] adding crossref type: %r',
                     (directivename, rolename, indextemplate, ref_nodeclass, objname))

        # create a subclass of Target as the new directive
        directive = type(directivename,  # type: ignore
                         (Target, object),
                         {'indextemplate': indextemplate})

        self.add_directive_to_domain('std', directivename, directive)
        self.add_role_to_domain('std', rolename, XRefRole(innernodeclass=ref_nodeclass))

        object_types = self.domain_object_types.setdefault('std', {})
        object_types[directivename] = ObjType(objname or directivename, rolename)

    def add_source_parser(self, suffix, parser):
        # type: (unicode, Type[Parser]) -> None
        logger.debug('[app] adding search source_parser: %r, %r', suffix, parser)
        if suffix in self.source_parsers:
            raise ExtensionError(__('source_parser for %r is already registered') % suffix)
        self.source_parsers[suffix] = parser

    def get_source_parser(self, filename):
        # type: (unicode) -> Type[Parser]
        for suffix, parser_class in iteritems(self.source_parsers):
            if filename.endswith(suffix):
                break
        else:
            # use special parser for unknown file-extension '*' (if exists)
            parser_class = self.source_parsers.get('*')

        if parser_class is None:
            raise SphinxError(__('source_parser for %s not registered') % filename)
        else:
            if isinstance(parser_class, string_types):
                parser_class = import_object(parser_class, 'source parser')  # type: ignore
            return parser_class

    def get_source_parsers(self):
        # type: () -> Dict[unicode, Parser]
        return self.source_parsers

    def create_source_parser(self, app, filename):
        # type: (Sphinx, unicode) -> Parser
        parser_class = self.get_source_parser(filename)
        parser = parser_class()
        if isinstance(parser, SphinxParser):
            parser.set_application(app)
        return parser

    def add_source_input(self, input_class):
        # type: (Type[Input]) -> None
        for filetype in input_class.supported:
            if filetype in self.source_inputs:
                raise ExtensionError(__('source_input for %r is already registered') %
                                     filetype)
            self.source_inputs[filetype] = input_class

    def get_source_input(self, filename):
        # type: (unicode) -> Type[Input]
        parser = self.get_source_parser(filename)
        for filetype in parser.supported:
            if filetype in self.source_inputs:
                input_class = self.source_inputs[filetype]
                break
        else:
            # use special source_input for unknown file-type '*' (if exists)
            input_class = self.source_inputs.get('*')

        if input_class is None:
            raise SphinxError(__('source_input for %s not registered') % filename)
        else:
            return input_class

    def add_translator(self, name, translator):
        # type: (unicode, Type[nodes.NodeVisitor]) -> None
        logger.info(bold(__('Change of translator for the %s builder.') % name))
        self.translators[name] = translator

    def get_translator_class(self, builder):
        # type: (Builder) -> Type[nodes.NodeVisitor]
        return self.translators.get(builder.name,
                                    builder.default_translator_class)

    def create_translator(self, builder, document):
        # type: (Builder, nodes.Node) -> nodes.NodeVisitor
        translator_class = self.get_translator_class(builder)
        return translator_class(builder, document)

    def add_transform(self, transform):
        # type: (Type[Transform]) -> None
        logger.debug('[app] adding transform: %r', transform)
        self.transforms.append(transform)

    def get_transforms(self):
        # type: () -> List[Type[Transform]]
        return self.transforms

    def add_post_transform(self, transform):
        # type: (Type[Transform]) -> None
        logger.debug('[app] adding post transform: %r', transform)
        self.post_transforms.append(transform)

    def get_post_transforms(self):
        # type: () -> List[Type[Transform]]
        return self.post_transforms

    def add_documenter(self, objtype, documenter):
        # type: (unicode, Type[Documenter]) -> None
        self.documenters[objtype] = documenter

    def add_autodoc_attrgetter(self, typ, attrgetter):
        # type: (Type, Callable[[Any, unicode, Any], Any]) -> None
        self.autodoc_attrgettrs[typ] = attrgetter

    def load_extension(self, app, extname):
        # type: (Sphinx, unicode) -> None
        """Load a Sphinx extension."""
        if extname in app.extensions:  # alread loaded
            return
        if extname in EXTENSION_BLACKLIST:
            logger.warning(__('the extension %r was already merged with Sphinx since '
                              'version %s; this extension is ignored.'),
                           extname, EXTENSION_BLACKLIST[extname])
            return

        # update loading context
        app._setting_up_extension.append(extname)

        try:
            mod = __import__(extname, None, None, ['setup'])
        except ImportError as err:
            logger.verbose(__('Original exception:\n') + traceback.format_exc())
            raise ExtensionError(__('Could not import extension %s') % extname, err)

        if not hasattr(mod, 'setup'):
            logger.warning(__('extension %r has no setup() function; is it really '
                              'a Sphinx extension module?'), extname)
            metadata = {}  # type: Dict[unicode, Any]
        else:
            try:
                metadata = mod.setup(app)
            except VersionRequirementError as err:
                # add the extension name to the version required
                raise VersionRequirementError(
                    __('The %s extension used by this project needs at least '
                       'Sphinx v%s; it therefore cannot be built with this '
                       'version.') % (extname, err)
                )

        if metadata is None:
            metadata = {}
            if extname == 'rst2pdf.pdfbuilder':
                metadata['parallel_read_safe'] = True
        elif not isinstance(metadata, dict):
            logger.warning(__('extension %r returned an unsupported object from '
                              'its setup() function; it should return None or a '
                              'metadata dictionary'), extname)
            metadata = {}

        app.extensions[extname] = Extension(extname, mod, **metadata)
        app._setting_up_extension.pop()
