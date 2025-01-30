"""Insert links to objects documented in remote Sphinx documentation.

This works as follows:

* Each Sphinx HTML build creates a file named "objects.inv" that contains a
  mapping from object names to URIs relative to the HTML set's root.

* Projects using the Intersphinx extension can specify links to such mapping
  files in the `intersphinx_mapping` config value.  The mapping will then be
  used to resolve otherwise missing references to objects into links to the
  other documentation.

* By default, the mapping file is assumed to be at the same location as the
  rest of the documentation; however, the location of the mapping file can
  also be specified individually, e.g. if the docs should be buildable
  without Internet access.
"""

from __future__ import annotations

__all__ = (
    'InventoryAdapter',
    'fetch_inventory',
    'load_mappings',
    'validate_intersphinx_mapping',
    'IntersphinxRoleResolver',
    'inventory_exists',
    'install_dispatcher',
    'resolve_reference_in_inventory',
    'resolve_reference_any_inventory',
    'resolve_reference_detect_inventory',
    'missing_reference',
    'IntersphinxDispatcher',
    'IntersphinxRole',
    'inspect_main',
)

from typing import TYPE_CHECKING

import sphinx
from sphinx.ext.intersphinx._cli import inspect_main
from sphinx.ext.intersphinx._load import (
    fetch_inventory,
    load_mappings,
    validate_intersphinx_mapping,
)
from sphinx.ext.intersphinx._resolve import (
    IntersphinxDispatcher,
    IntersphinxRole,
    IntersphinxRoleResolver,
    install_dispatcher,
    inventory_exists,
    missing_reference,
    resolve_reference_any_inventory,
    resolve_reference_detect_inventory,
    resolve_reference_in_inventory,
)
from sphinx.ext.intersphinx._shared import InventoryAdapter

if TYPE_CHECKING:
    from sphinx.application import Sphinx
    from sphinx.util.typing import ExtensionMetadata


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_config_value('intersphinx_mapping', {}, 'env')
    app.add_config_value('intersphinx_cache_limit', 5, '')
    app.add_config_value('intersphinx_timeout', None, '')
    app.add_config_value('intersphinx_disabled_reftypes', ['std:doc'], 'env')
    app.connect('config-inited', validate_intersphinx_mapping, priority=800)
    app.connect('builder-inited', load_mappings)
    app.connect('source-read', install_dispatcher)
    app.connect('missing-reference', missing_reference)
    app.add_post_transform(IntersphinxRoleResolver)
    return {
        'version': sphinx.__display_version__,
        'env_version': 1,
        'parallel_read_safe': True,
    }


# deprecated name -> (object to return, canonical path or empty string, removal version)
_DEPRECATED_OBJECTS: dict[str, tuple[object, str, tuple[int, int]]] = {
    'normalize_intersphinx_mapping': (
        validate_intersphinx_mapping,
        'sphinx.ext.intersphinx.validate_intersphinx_mapping',
        (10, 0),
    ),
}


def __getattr__(name: str) -> object:
    if name not in _DEPRECATED_OBJECTS:
        msg = f'module {__name__!r} has no attribute {name!r}'
        raise AttributeError(msg)

    from sphinx.deprecation import _deprecation_warning

    deprecated_object, canonical_name, remove = _DEPRECATED_OBJECTS[name]
    _deprecation_warning(__name__, name, canonical_name, remove=remove)
    return deprecated_object
