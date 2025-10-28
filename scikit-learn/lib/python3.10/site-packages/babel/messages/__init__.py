"""
    babel.messages
    ~~~~~~~~~~~~~~

    Support for ``gettext`` message catalogs.

    :copyright: (c) 2013-2025 by the Babel Team.
    :license: BSD, see LICENSE for more details.
"""

from babel.messages.catalog import (
    Catalog,
    Message,
    TranslationError,
)

__all__ = [
    "Catalog",
    "Message",
    "TranslationError",
]
