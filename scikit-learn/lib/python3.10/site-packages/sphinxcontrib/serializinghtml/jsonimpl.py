"""JSON serializer implementation wrapper."""

from __future__ import annotations

import json
from collections import UserString
from typing import IO, Any


class SphinxJSONEncoder(json.JSONEncoder):
    """JSONEncoder subclass that forces translation proxies."""
    def default(self, obj: Any) -> str:
        if isinstance(obj, UserString):
            return str(obj)
        return super().default(obj)


def dump(obj: Any, file: IO[str] | IO[bytes], *args: Any, **kwds: Any) -> None:
    kwds['cls'] = SphinxJSONEncoder
    json.dump(obj, file, *args, **kwds)


def dumps(obj: Any, *args: Any, **kwds: Any) -> str:
    kwds['cls'] = SphinxJSONEncoder
    return json.dumps(obj, *args, **kwds)


def load(*args: Any, **kwds: Any) -> Any:
    return json.load(*args, **kwds)


def loads(*args: Any, **kwds: Any) -> Any:
    return json.loads(*args, **kwds)
