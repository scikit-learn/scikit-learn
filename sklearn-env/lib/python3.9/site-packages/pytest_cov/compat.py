try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pytest

StringIO  # pyflakes, this is for re-export


if hasattr(pytest, 'hookimpl'):
    hookwrapper = pytest.hookimpl(hookwrapper=True)
else:
    hookwrapper = pytest.mark.hookwrapper


class SessionWrapper:
    def __init__(self, session):
        self._session = session
        if hasattr(session, 'testsfailed'):
            self._attr = 'testsfailed'
        else:
            self._attr = '_testsfailed'

    @property
    def testsfailed(self):
        return getattr(self._session, self._attr)

    @testsfailed.setter
    def testsfailed(self, value):
        setattr(self._session, self._attr, value)
