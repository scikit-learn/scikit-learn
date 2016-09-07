# -*- coding: utf-8 -*-
"""
    linkcode
    ~~~~~~~~

    Add external links to module code in Python object descriptions.

    :copyright: Copyright 2007-2011 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.

"""
from __future__ import division, absolute_import, print_function

import warnings
import collections

warnings.warn("This extension has been accepted to Sphinx upstream. "
              "Use the version from there (Sphinx >= 1.2) "
              "https://bitbucket.org/birkenfeld/sphinx/pull-request/47/sphinxextlinkcode",
              FutureWarning, stacklevel=1)


from docutils import nodes

from sphinx import addnodes
from sphinx.locale import _
from sphinx.errors import SphinxError

class LinkcodeError(SphinxError):
    category = "linkcode error"

def doctree_read(app, doctree):
    env = app.builder.env

    resolve_target = getattr(env.config, 'linkcode_resolve', None)
    if not isinstance(env.config.linkcode_resolve, collections.Callable):
        raise LinkcodeError(
            "Function `linkcode_resolve` is not given in conf.py")

    domain_keys = dict(
        py=['module', 'fullname'],
        c=['names'],
        cpp=['names'],
        js=['object', 'fullname'],
    )

    for objnode in doctree.traverse(addnodes.desc):
        domain = objnode.get('domain')
        uris = set()
        for signode in objnode:
            if not isinstance(signode, addnodes.desc_signature):
                continue

            # Convert signode to a specified format
            info = {}
            for key in domain_keys.get(domain, []):
                value = signode.get(key)
                if not value:
                    value = ''
                info[key] = value
            if not info:
                continue

            # Call user code to resolve the link
            uri = resolve_target(domain, info)
            if not uri:
                # no source
                continue

            if uri in uris or not uri:
                # only one link per name, please
                continue
            uris.add(uri)

            onlynode = addnodes.only(expr='html')
            onlynode += nodes.reference('', '', internal=False, refuri=uri)
            onlynode[0] += nodes.inline('', _('[source]'),
                                        classes=['viewcode-link'])
            signode += onlynode

def setup(app):
    app.connect('doctree-read', doctree_read)
    app.add_config_value('linkcode_resolve', None, '')
