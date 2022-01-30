"""
    sphinx.util.jsonimpl
    ~~~~~~~~~~~~~~~~~~~~

    JSON serializer implementation wrapper.

    :copyright: Copyright 2007-2019 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import json
from collections import UserString

if False:
    # For type annotation
    from typing import Any, IO  # NOQA


class SphinxJSONEncoder(json.JSONEncoder):
    """JSONEncoder subclass that forces translation proxies."""
    def default(self, obj):
        # type: (Any) -> str
        if isinstance(obj, UserString):
            return str(obj)
        return super().default(obj)


def dump(obj, fp, *args, **kwds):
    # type: (Any, IO, Any, Any) -> None
    kwds['cls'] = SphinxJSONEncoder
    json.dump(obj, fp, *args, **kwds)


def dumps(obj, *args, **kwds):
    # type: (Any, Any, Any) -> str
    kwds['cls'] = SphinxJSONEncoder
    return json.dumps(obj, *args, **kwds)


def load(*args, **kwds):
    # type: (Any, Any) -> Any
    return json.load(*args, **kwds)


def loads(*args, **kwds):
    # type: (Any, Any) -> Any
    return json.loads(*args, **kwds)
