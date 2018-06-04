from __future__ import absolute_import, division, print_function
import os

import six
import py
import attr

import _pytest
import _pytest._code

from _pytest.mark.structures import NodeKeywords

SEP = "/"

tracebackcutdir = py.path.local(_pytest.__file__).dirpath()


def _splitnode(nodeid):
    """Split a nodeid into constituent 'parts'.

    Node IDs are strings, and can be things like:
        ''
        'testing/code'
        'testing/code/test_excinfo.py'
        'testing/code/test_excinfo.py::TestFormattedExcinfo::()'

    Return values are lists e.g.
        []
        ['testing', 'code']
        ['testing', 'code', 'test_excinfo.py']
        ['testing', 'code', 'test_excinfo.py', 'TestFormattedExcinfo', '()']
    """
    if nodeid == '':
        # If there is no root node at all, return an empty list so the caller's logic can remain sane
        return []
    parts = nodeid.split(SEP)
    # Replace single last element 'test_foo.py::Bar::()' with multiple elements 'test_foo.py', 'Bar', '()'
    parts[-1:] = parts[-1].split("::")
    return parts


def ischildnode(baseid, nodeid):
    """Return True if the nodeid is a child node of the baseid.

    E.g. 'foo/bar::Baz::()' is a child of 'foo', 'foo/bar' and 'foo/bar::Baz', but not of 'foo/blorp'
    """
    base_parts = _splitnode(baseid)
    node_parts = _splitnode(nodeid)
    if len(node_parts) < len(base_parts):
        return False
    return node_parts[:len(base_parts)] == base_parts


@attr.s
class _CompatProperty(object):
    name = attr.ib()

    def __get__(self, obj, owner):
        if obj is None:
            return self

        # TODO: reenable in the features branch
        # warnings.warn(
        #     "usage of {owner!r}.{name} is deprecated, please use pytest.{name} instead".format(
        #         name=self.name, owner=type(owner).__name__),
        #     PendingDeprecationWarning, stacklevel=2)
        return getattr(__import__('pytest'), self.name)


class Node(object):
    """ base class for Collector and Item the test collection tree.
    Collector subclasses have children, Items are terminal nodes."""

    def __init__(self, name, parent=None, config=None, session=None, fspath=None, nodeid=None):
        #: a unique name within the scope of the parent node
        self.name = name

        #: the parent collector node.
        self.parent = parent

        #: the pytest config object
        self.config = config or parent.config

        #: the session this node is part of
        self.session = session or parent.session

        #: filesystem path where this node was collected from (can be None)
        self.fspath = fspath or getattr(parent, 'fspath', None)

        #: keywords/markers collected from all scopes
        self.keywords = NodeKeywords(self)

        #: allow adding of extra keywords to use for matching
        self.extra_keyword_matches = set()

        # used for storing artificial fixturedefs for direct parametrization
        self._name2pseudofixturedef = {}

        if nodeid is not None:
            self._nodeid = nodeid
        else:
            assert parent is not None
            self._nodeid = self.parent.nodeid + "::" + self.name

    @property
    def ihook(self):
        """ fspath sensitive hook proxy used to call pytest hooks"""
        return self.session.gethookproxy(self.fspath)

    Module = _CompatProperty("Module")
    Class = _CompatProperty("Class")
    Instance = _CompatProperty("Instance")
    Function = _CompatProperty("Function")
    File = _CompatProperty("File")
    Item = _CompatProperty("Item")

    def _getcustomclass(self, name):
        maybe_compatprop = getattr(type(self), name)
        if isinstance(maybe_compatprop, _CompatProperty):
            return getattr(__import__('pytest'), name)
        else:
            cls = getattr(self, name)
            # TODO: reenable in the features branch
            # warnings.warn("use of node.%s is deprecated, "
            #    "use pytest_pycollect_makeitem(...) to create custom "
            #    "collection nodes" % name, category=DeprecationWarning)
        return cls

    def __repr__(self):
        return "<%s %r>" % (self.__class__.__name__,
                            getattr(self, 'name', None))

    def warn(self, code, message):
        """ generate a warning with the given code and message for this
        item. """
        assert isinstance(code, str)
        fslocation = getattr(self, "location", None)
        if fslocation is None:
            fslocation = getattr(self, "fspath", None)
        self.ihook.pytest_logwarning.call_historic(kwargs=dict(
            code=code, message=message,
            nodeid=self.nodeid, fslocation=fslocation))

    # methods for ordering nodes
    @property
    def nodeid(self):
        """ a ::-separated string denoting its collection tree address. """
        return self._nodeid

    def __hash__(self):
        return hash(self.nodeid)

    def setup(self):
        pass

    def teardown(self):
        pass

    def listchain(self):
        """ return list of all parent collectors up to self,
            starting from root of collection tree. """
        chain = []
        item = self
        while item is not None:
            chain.append(item)
            item = item.parent
        chain.reverse()
        return chain

    def add_marker(self, marker):
        """ dynamically add a marker object to the node.

        ``marker`` can be a string or pytest.mark.* instance.
        """
        from _pytest.mark import MarkDecorator, MARK_GEN
        if isinstance(marker, six.string_types):
            marker = getattr(MARK_GEN, marker)
        elif not isinstance(marker, MarkDecorator):
            raise ValueError("is not a string or pytest.mark.* Marker")
        self.keywords[marker.name] = marker

    def get_marker(self, name):
        """ get a marker object from this node or None if
        the node doesn't have a marker with that name. """
        val = self.keywords.get(name, None)
        if val is not None:
            from _pytest.mark import MarkInfo, MarkDecorator
            if isinstance(val, (MarkDecorator, MarkInfo)):
                return val

    def listextrakeywords(self):
        """ Return a set of all extra keywords in self and any parents."""
        extra_keywords = set()
        for item in self.listchain():
            extra_keywords.update(item.extra_keyword_matches)
        return extra_keywords

    def listnames(self):
        return [x.name for x in self.listchain()]

    def addfinalizer(self, fin):
        """ register a function to be called when this node is finalized.

        This method can only be called when this node is active
        in a setup chain, for example during self.setup().
        """
        self.session._setupstate.addfinalizer(fin, self)

    def getparent(self, cls):
        """ get the next parent node (including ourself)
        which is an instance of the given class"""
        current = self
        while current and not isinstance(current, cls):
            current = current.parent
        return current

    def _prunetraceback(self, excinfo):
        pass

    def _repr_failure_py(self, excinfo, style=None):
        fm = self.session._fixturemanager
        if excinfo.errisinstance(fm.FixtureLookupError):
            return excinfo.value.formatrepr()
        tbfilter = True
        if self.config.option.fulltrace:
            style = "long"
        else:
            tb = _pytest._code.Traceback([excinfo.traceback[-1]])
            self._prunetraceback(excinfo)
            if len(excinfo.traceback) == 0:
                excinfo.traceback = tb
            tbfilter = False  # prunetraceback already does it
            if style == "auto":
                style = "long"
        # XXX should excinfo.getrepr record all data and toterminal() process it?
        if style is None:
            if self.config.option.tbstyle == "short":
                style = "short"
            else:
                style = "long"

        try:
            os.getcwd()
            abspath = False
        except OSError:
            abspath = True

        return excinfo.getrepr(funcargs=True, abspath=abspath,
                               showlocals=self.config.option.showlocals,
                               style=style, tbfilter=tbfilter)

    repr_failure = _repr_failure_py


class Collector(Node):
    """ Collector instances create children through collect()
        and thus iteratively build a tree.
    """

    class CollectError(Exception):
        """ an error during collection, contains a custom message. """

    def collect(self):
        """ returns a list of children (items and collectors)
            for this collection node.
        """
        raise NotImplementedError("abstract")

    def repr_failure(self, excinfo):
        """ represent a collection failure. """
        if excinfo.errisinstance(self.CollectError):
            exc = excinfo.value
            return str(exc.args[0])
        return self._repr_failure_py(excinfo, style="short")

    def _prunetraceback(self, excinfo):
        if hasattr(self, 'fspath'):
            traceback = excinfo.traceback
            ntraceback = traceback.cut(path=self.fspath)
            if ntraceback == traceback:
                ntraceback = ntraceback.cut(excludepath=tracebackcutdir)
            excinfo.traceback = ntraceback.filter()


def _check_initialpaths_for_relpath(session, fspath):
    for initial_path in session._initialpaths:
        if fspath.common(initial_path) == initial_path:
            return fspath.relto(initial_path.dirname)


class FSCollector(Collector):
    def __init__(self, fspath, parent=None, config=None, session=None, nodeid=None):
        fspath = py.path.local(fspath)  # xxx only for test_resultlog.py?
        name = fspath.basename
        if parent is not None:
            rel = fspath.relto(parent.fspath)
            if rel:
                name = rel
            name = name.replace(os.sep, SEP)
        self.fspath = fspath

        session = session or parent.session

        if nodeid is None:
            nodeid = self.fspath.relto(session.config.rootdir)

            if not nodeid:
                nodeid = _check_initialpaths_for_relpath(session, fspath)
            if os.sep != SEP:
                nodeid = nodeid.replace(os.sep, SEP)

        super(FSCollector, self).__init__(name, parent, config, session, nodeid=nodeid, fspath=fspath)


class File(FSCollector):
    """ base class for collecting tests from a file. """


class Item(Node):
    """ a basic test invocation item. Note that for a single function
    there might be multiple test invocation items.
    """
    nextitem = None

    def __init__(self, name, parent=None, config=None, session=None, nodeid=None):
        super(Item, self).__init__(name, parent, config, session, nodeid=nodeid)
        self._report_sections = []

        #: user properties is a list of tuples (name, value) that holds user
        #: defined properties for this test.
        self.user_properties = []

    def add_report_section(self, when, key, content):
        """
        Adds a new report section, similar to what's done internally to add stdout and
        stderr captured output::

            item.add_report_section("call", "stdout", "report section contents")

        :param str when:
            One of the possible capture states, ``"setup"``, ``"call"``, ``"teardown"``.
        :param str key:
            Name of the section, can be customized at will. Pytest uses ``"stdout"`` and
            ``"stderr"`` internally.

        :param str content:
            The full contents as a string.
        """
        if content:
            self._report_sections.append((when, key, content))

    def reportinfo(self):
        return self.fspath, None, ""

    @property
    def location(self):
        try:
            return self._location
        except AttributeError:
            location = self.reportinfo()
            # bestrelpath is a quite slow function
            cache = self.config.__dict__.setdefault("_bestrelpathcache", {})
            try:
                fspath = cache[location[0]]
            except KeyError:
                fspath = self.session.fspath.bestrelpath(location[0])
                cache[location[0]] = fspath
            location = (fspath, location[1], str(location[2]))
            self._location = location
            return location
