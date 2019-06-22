# coding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

import datetime
import copy

# ToDo: at least on PY3 you could probably attach the tzinfo correctly to the object
#       a more complete datetime might be used by safe loading as well

if False:  # MYPY
    from typing import Any, Dict, Optional, List  # NOQA


class TimeStamp(datetime.datetime):
    def __init__(self, *args, **kw):
        # type: (Any, Any) -> None
        self._yaml = dict(t=False, tz=None, delta=0)  # type: Dict[Any, Any]

    def __new__(cls, *args, **kw):  # datetime is immutable
        # type: (Any, Any) -> Any
        return datetime.datetime.__new__(cls, *args, **kw)  # type: ignore

    def __deepcopy__(self, memo):
        # type: (Any) -> Any
        ts = TimeStamp(self.year, self.month, self.day, self.hour, self.minute, self.second)
        ts._yaml = copy.deepcopy(self._yaml)
        return ts
