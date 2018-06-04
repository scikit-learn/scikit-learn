""" support for providing temporary directories to test functions.  """
from __future__ import absolute_import, division, print_function

import re

import pytest
import py
from _pytest.monkeypatch import MonkeyPatch


class TempdirFactory(object):
    """Factory for temporary directories under the common base temp directory.

    The base directory can be configured using the ``--basetemp`` option.
    """

    def __init__(self, config):
        self.config = config
        self.trace = config.trace.get("tmpdir")

    def ensuretemp(self, string, dir=1):
        """ (deprecated) return temporary directory path with
            the given string as the trailing part.  It is usually
            better to use the 'tmpdir' function argument which
            provides an empty unique-per-test-invocation directory
            and is guaranteed to be empty.
        """
        # py.log._apiwarn(">1.1", "use tmpdir function argument")
        return self.getbasetemp().ensure(string, dir=dir)

    def mktemp(self, basename, numbered=True):
        """Create a subdirectory of the base temporary directory and return it.
        If ``numbered``, ensure the directory is unique by adding a number
        prefix greater than any existing one.
        """
        basetemp = self.getbasetemp()
        if not numbered:
            p = basetemp.mkdir(basename)
        else:
            p = py.path.local.make_numbered_dir(prefix=basename,
                                                keep=0, rootdir=basetemp, lock_timeout=None)
        self.trace("mktemp", p)
        return p

    def getbasetemp(self):
        """ return base temporary directory. """
        try:
            return self._basetemp
        except AttributeError:
            basetemp = self.config.option.basetemp
            if basetemp:
                basetemp = py.path.local(basetemp)
                if basetemp.check():
                    basetemp.remove()
                basetemp.mkdir()
            else:
                temproot = py.path.local.get_temproot()
                user = get_user()
                if user:
                    # use a sub-directory in the temproot to speed-up
                    # make_numbered_dir() call
                    rootdir = temproot.join('pytest-of-%s' % user)
                else:
                    rootdir = temproot
                rootdir.ensure(dir=1)
                basetemp = py.path.local.make_numbered_dir(prefix='pytest-',
                                                           rootdir=rootdir)
            self._basetemp = t = basetemp.realpath()
            self.trace("new basetemp", t)
            return t

    def finish(self):
        self.trace("finish")


def get_user():
    """Return the current user name, or None if getuser() does not work
    in the current environment (see #1010).
    """
    import getpass
    try:
        return getpass.getuser()
    except (ImportError, KeyError):
        return None


# backward compatibility
TempdirHandler = TempdirFactory


def pytest_configure(config):
    """Create a TempdirFactory and attach it to the config object.

    This is to comply with existing plugins which expect the handler to be
    available at pytest_configure time, but ideally should be moved entirely
    to the tmpdir_factory session fixture.
    """
    mp = MonkeyPatch()
    t = TempdirFactory(config)
    config._cleanup.extend([mp.undo, t.finish])
    mp.setattr(config, '_tmpdirhandler', t, raising=False)
    mp.setattr(pytest, 'ensuretemp', t.ensuretemp, raising=False)


@pytest.fixture(scope='session')
def tmpdir_factory(request):
    """Return a TempdirFactory instance for the test session.
    """
    return request.config._tmpdirhandler


@pytest.fixture
def tmpdir(request, tmpdir_factory):
    """Return a temporary directory path object
    which is unique to each test function invocation,
    created as a sub directory of the base temporary
    directory.  The returned object is a `py.path.local`_
    path object.

    .. _`py.path.local`: https://py.readthedocs.io/en/latest/path.html
    """
    name = request.node.name
    name = re.sub(r"[\W]", "_", name)
    MAXVAL = 30
    if len(name) > MAXVAL:
        name = name[:MAXVAL]
    x = tmpdir_factory.mktemp(name, numbered=True)
    return x
