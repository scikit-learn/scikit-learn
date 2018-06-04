from __future__ import absolute_import, division, print_function

import functools
import inspect
import sys
import warnings
from collections import OrderedDict, deque, defaultdict

import attr
import py
from py._code.code import FormattedExcinfo

import _pytest
from _pytest import nodes
from _pytest._code.code import TerminalRepr
from _pytest.compat import (
    NOTSET, exc_clear, _format_args,
    getfslineno, get_real_func,
    is_generator, isclass, getimfunc,
    getlocation, getfuncargnames,
    safe_getattr,
    FuncargnamesCompatAttr,
)
from _pytest.outcomes import fail, TEST_OUTCOME


@attr.s(frozen=True)
class PseudoFixtureDef(object):
    cached_result = attr.ib()
    scope = attr.ib()


def pytest_sessionstart(session):
    import _pytest.python
    import _pytest.nodes

    scopename2class.update({
        'class': _pytest.python.Class,
        'module': _pytest.python.Module,
        'function': _pytest.nodes.Item,
        'session': _pytest.main.Session,
    })
    session._fixturemanager = FixtureManager(session)


scopename2class = {}


scope2props = dict(session=())
scope2props["module"] = ("fspath", "module")
scope2props["class"] = scope2props["module"] + ("cls",)
scope2props["instance"] = scope2props["class"] + ("instance", )
scope2props["function"] = scope2props["instance"] + ("function", "keywords")


def scopeproperty(name=None, doc=None):
    def decoratescope(func):
        scopename = name or func.__name__

        def provide(self):
            if func.__name__ in scope2props[self.scope]:
                return func(self)
            raise AttributeError("%s not available in %s-scoped context" % (
                scopename, self.scope))

        return property(provide, None, None, func.__doc__)
    return decoratescope


def get_scope_node(node, scope):
    cls = scopename2class.get(scope)
    if cls is None:
        raise ValueError("unknown scope")
    return node.getparent(cls)


def add_funcarg_pseudo_fixture_def(collector, metafunc, fixturemanager):
    # this function will transform all collected calls to a functions
    # if they use direct funcargs (i.e. direct parametrization)
    # because we want later test execution to be able to rely on
    # an existing FixtureDef structure for all arguments.
    # XXX we can probably avoid this algorithm  if we modify CallSpec2
    # to directly care for creating the fixturedefs within its methods.
    if not metafunc._calls[0].funcargs:
        return  # this function call does not have direct parametrization
    # collect funcargs of all callspecs into a list of values
    arg2params = {}
    arg2scope = {}
    for callspec in metafunc._calls:
        for argname, argvalue in callspec.funcargs.items():
            assert argname not in callspec.params
            callspec.params[argname] = argvalue
            arg2params_list = arg2params.setdefault(argname, [])
            callspec.indices[argname] = len(arg2params_list)
            arg2params_list.append(argvalue)
            if argname not in arg2scope:
                scopenum = callspec._arg2scopenum.get(argname,
                                                      scopenum_function)
                arg2scope[argname] = scopes[scopenum]
        callspec.funcargs.clear()

    # register artificial FixtureDef's so that later at test execution
    # time we can rely on a proper FixtureDef to exist for fixture setup.
    arg2fixturedefs = metafunc._arg2fixturedefs
    for argname, valuelist in arg2params.items():
        # if we have a scope that is higher than function we need
        # to make sure we only ever create an according fixturedef on
        # a per-scope basis. We thus store and cache the fixturedef on the
        # node related to the scope.
        scope = arg2scope[argname]
        node = None
        if scope != "function":
            node = get_scope_node(collector, scope)
            if node is None:
                assert scope == "class" and isinstance(collector, _pytest.python.Module)
                # use module-level collector for class-scope (for now)
                node = collector
        if node and argname in node._name2pseudofixturedef:
            arg2fixturedefs[argname] = [node._name2pseudofixturedef[argname]]
        else:
            fixturedef = FixtureDef(fixturemanager, '', argname,
                                    get_direct_param_fixture_func,
                                    arg2scope[argname],
                                    valuelist, False, False)
            arg2fixturedefs[argname] = [fixturedef]
            if node is not None:
                node._name2pseudofixturedef[argname] = fixturedef


def getfixturemarker(obj):
    """ return fixturemarker or None if it doesn't exist or raised
    exceptions."""
    try:
        return getattr(obj, "_pytestfixturefunction", None)
    except TEST_OUTCOME:
        # some objects raise errors like request (from flask import request)
        # we don't expect them to be fixture functions
        return None


def get_parametrized_fixture_keys(item, scopenum):
    """ return list of keys for all parametrized arguments which match
    the specified scope. """
    assert scopenum < scopenum_function  # function
    try:
        cs = item.callspec
    except AttributeError:
        pass
    else:
        # cs.indices.items() is random order of argnames.  Need to
        # sort this so that different calls to
        # get_parametrized_fixture_keys will be deterministic.
        for argname, param_index in sorted(cs.indices.items()):
            if cs._arg2scopenum[argname] != scopenum:
                continue
            if scopenum == 0:    # session
                key = (argname, param_index)
            elif scopenum == 1:  # module
                key = (argname, param_index, item.fspath)
            elif scopenum == 2:  # class
                key = (argname, param_index, item.fspath, item.cls)
            yield key


# algorithm for sorting on a per-parametrized resource setup basis
# it is called for scopenum==0 (session) first and performs sorting
# down to the lower scopes such as to minimize number of "high scope"
# setups and teardowns

def reorder_items(items):
    argkeys_cache = {}
    items_by_argkey = {}
    for scopenum in range(0, scopenum_function):
        argkeys_cache[scopenum] = d = {}
        items_by_argkey[scopenum] = item_d = defaultdict(deque)
        for item in items:
            keys = OrderedDict.fromkeys(get_parametrized_fixture_keys(item, scopenum))
            if keys:
                d[item] = keys
                for key in keys:
                    item_d[key].append(item)
    items = OrderedDict.fromkeys(items)
    return list(reorder_items_atscope(items, argkeys_cache, items_by_argkey, 0))


def fix_cache_order(item, argkeys_cache, items_by_argkey):
    for scopenum in range(0, scopenum_function):
        for key in argkeys_cache[scopenum].get(item, []):
            items_by_argkey[scopenum][key].appendleft(item)


def reorder_items_atscope(items, argkeys_cache, items_by_argkey, scopenum):
    if scopenum >= scopenum_function or len(items) < 3:
        return items
    ignore = set()
    items_deque = deque(items)
    items_done = OrderedDict()
    scoped_items_by_argkey = items_by_argkey[scopenum]
    scoped_argkeys_cache = argkeys_cache[scopenum]
    while items_deque:
        no_argkey_group = OrderedDict()
        slicing_argkey = None
        while items_deque:
            item = items_deque.popleft()
            if item in items_done or item in no_argkey_group:
                continue
            argkeys = OrderedDict.fromkeys(k for k in scoped_argkeys_cache.get(item, []) if k not in ignore)
            if not argkeys:
                no_argkey_group[item] = None
            else:
                slicing_argkey, _ = argkeys.popitem()
                # we don't have to remove relevant items from later in the deque because they'll just be ignored
                matching_items = [i for i in scoped_items_by_argkey[slicing_argkey] if i in items]
                for i in reversed(matching_items):
                    fix_cache_order(i, argkeys_cache, items_by_argkey)
                    items_deque.appendleft(i)
                break
        if no_argkey_group:
            no_argkey_group = reorder_items_atscope(
                                no_argkey_group, argkeys_cache, items_by_argkey, scopenum + 1)
            for item in no_argkey_group:
                items_done[item] = None
        ignore.add(slicing_argkey)
    return items_done


def fillfixtures(function):
    """ fill missing funcargs for a test function. """
    try:
        request = function._request
    except AttributeError:
        # XXX this special code path is only expected to execute
        # with the oejskit plugin.  It uses classes with funcargs
        # and we thus have to work a bit to allow this.
        fm = function.session._fixturemanager
        fi = fm.getfixtureinfo(function.parent, function.obj, None)
        function._fixtureinfo = fi
        request = function._request = FixtureRequest(function)
        request._fillfixtures()
        # prune out funcargs for jstests
        newfuncargs = {}
        for name in fi.argnames:
            newfuncargs[name] = function.funcargs[name]
        function.funcargs = newfuncargs
    else:
        request._fillfixtures()


def get_direct_param_fixture_func(request):
    return request.param


class FuncFixtureInfo(object):
    def __init__(self, argnames, names_closure, name2fixturedefs):
        self.argnames = argnames
        self.names_closure = names_closure
        self.name2fixturedefs = name2fixturedefs


class FixtureRequest(FuncargnamesCompatAttr):
    """ A request for a fixture from a test or fixture function.

    A request object gives access to the requesting test context
    and has an optional ``param`` attribute in case
    the fixture is parametrized indirectly.
    """

    def __init__(self, pyfuncitem):
        self._pyfuncitem = pyfuncitem
        #: fixture for which this request is being performed
        self.fixturename = None
        #: Scope string, one of "function", "class", "module", "session"
        self.scope = "function"
        self._fixture_defs = {}  # argname -> FixtureDef
        fixtureinfo = pyfuncitem._fixtureinfo
        self._arg2fixturedefs = fixtureinfo.name2fixturedefs.copy()
        self._arg2index = {}
        self._fixturemanager = pyfuncitem.session._fixturemanager

    @property
    def fixturenames(self):
        # backward incompatible note: now a readonly property
        return list(self._pyfuncitem._fixtureinfo.names_closure)

    @property
    def node(self):
        """ underlying collection node (depends on current request scope)"""
        return self._getscopeitem(self.scope)

    def _getnextfixturedef(self, argname):
        fixturedefs = self._arg2fixturedefs.get(argname, None)
        if fixturedefs is None:
            # we arrive here because of a  a dynamic call to
            # getfixturevalue(argname) usage which was naturally
            # not known at parsing/collection time
            parentid = self._pyfuncitem.parent.nodeid
            fixturedefs = self._fixturemanager.getfixturedefs(argname, parentid)
            self._arg2fixturedefs[argname] = fixturedefs
        # fixturedefs list is immutable so we maintain a decreasing index
        index = self._arg2index.get(argname, 0) - 1
        if fixturedefs is None or (-index > len(fixturedefs)):
            raise FixtureLookupError(argname, self)
        self._arg2index[argname] = index
        return fixturedefs[index]

    @property
    def config(self):
        """ the pytest config object associated with this request. """
        return self._pyfuncitem.config

    @scopeproperty()
    def function(self):
        """ test function object if the request has a per-function scope. """
        return self._pyfuncitem.obj

    @scopeproperty("class")
    def cls(self):
        """ class (can be None) where the test function was collected. """
        clscol = self._pyfuncitem.getparent(_pytest.python.Class)
        if clscol:
            return clscol.obj

    @property
    def instance(self):
        """ instance (can be None) on which test function was collected. """
        # unittest support hack, see _pytest.unittest.TestCaseFunction
        try:
            return self._pyfuncitem._testcase
        except AttributeError:
            function = getattr(self, "function", None)
            if function is not None:
                return py.builtin._getimself(function)

    @scopeproperty()
    def module(self):
        """ python module object where the test function was collected. """
        return self._pyfuncitem.getparent(_pytest.python.Module).obj

    @scopeproperty()
    def fspath(self):
        """ the file system path of the test module which collected this test. """
        return self._pyfuncitem.fspath

    @property
    def keywords(self):
        """ keywords/markers dictionary for the underlying node. """
        return self.node.keywords

    @property
    def session(self):
        """ pytest session object. """
        return self._pyfuncitem.session

    def addfinalizer(self, finalizer):
        """ add finalizer/teardown function to be called after the
        last test within the requesting test context finished
        execution. """
        # XXX usually this method is shadowed by fixturedef specific ones
        self._addfinalizer(finalizer, scope=self.scope)

    def _addfinalizer(self, finalizer, scope):
        colitem = self._getscopeitem(scope)
        self._pyfuncitem.session._setupstate.addfinalizer(
            finalizer=finalizer, colitem=colitem)

    def applymarker(self, marker):
        """ Apply a marker to a single test function invocation.
        This method is useful if you don't want to have a keyword/marker
        on all function invocations.

        :arg marker: a :py:class:`_pytest.mark.MarkDecorator` object
            created by a call to ``pytest.mark.NAME(...)``.
        """
        try:
            self.node.keywords[marker.markname] = marker
        except AttributeError:
            raise ValueError(marker)

    def raiseerror(self, msg):
        """ raise a FixtureLookupError with the given message. """
        raise self._fixturemanager.FixtureLookupError(None, self, msg)

    def _fillfixtures(self):
        item = self._pyfuncitem
        fixturenames = getattr(item, "fixturenames", self.fixturenames)
        for argname in fixturenames:
            if argname not in item.funcargs:
                item.funcargs[argname] = self.getfixturevalue(argname)

    def cached_setup(self, setup, teardown=None, scope="module", extrakey=None):
        """ (deprecated) Return a testing resource managed by ``setup`` &
        ``teardown`` calls.  ``scope`` and ``extrakey`` determine when the
        ``teardown`` function will be called so that subsequent calls to
        ``setup`` would recreate the resource.  With pytest-2.3 you often
        do not need ``cached_setup()`` as you can directly declare a scope
        on a fixture function and register a finalizer through
        ``request.addfinalizer()``.

        :arg teardown: function receiving a previously setup resource.
        :arg setup: a no-argument function creating a resource.
        :arg scope: a string value out of ``function``, ``class``, ``module``
            or ``session`` indicating the caching lifecycle of the resource.
        :arg extrakey: added to internal caching key of (funcargname, scope).
        """
        if not hasattr(self.config, '_setupcache'):
            self.config._setupcache = {}  # XXX weakref?
        cachekey = (self.fixturename, self._getscopeitem(scope), extrakey)
        cache = self.config._setupcache
        try:
            val = cache[cachekey]
        except KeyError:
            self._check_scope(self.fixturename, self.scope, scope)
            val = setup()
            cache[cachekey] = val
            if teardown is not None:
                def finalizer():
                    del cache[cachekey]
                    teardown(val)
                self._addfinalizer(finalizer, scope=scope)
        return val

    def getfixturevalue(self, argname):
        """ Dynamically run a named fixture function.

        Declaring fixtures via function argument is recommended where possible.
        But if you can only decide whether to use another fixture at test
        setup time, you may use this function to retrieve it inside a fixture
        or test function body.
        """
        return self._get_active_fixturedef(argname).cached_result[0]

    def getfuncargvalue(self, argname):
        """ Deprecated, use getfixturevalue. """
        from _pytest import deprecated
        warnings.warn(
            deprecated.GETFUNCARGVALUE,
            DeprecationWarning,
            stacklevel=2)
        return self.getfixturevalue(argname)

    def _get_active_fixturedef(self, argname):
        try:
            return self._fixture_defs[argname]
        except KeyError:
            try:
                fixturedef = self._getnextfixturedef(argname)
            except FixtureLookupError:
                if argname == "request":
                    cached_result = (self, [0], None)
                    scope = "function"
                    return PseudoFixtureDef(cached_result, scope)
                raise
        # remove indent to prevent the python3 exception
        # from leaking into the call
        self._compute_fixture_value(fixturedef)
        self._fixture_defs[argname] = fixturedef
        return fixturedef

    def _get_fixturestack(self):
        current = self
        values = []
        while 1:
            fixturedef = getattr(current, "_fixturedef", None)
            if fixturedef is None:
                values.reverse()
                return values
            values.append(fixturedef)
            current = current._parent_request

    def _compute_fixture_value(self, fixturedef):
        """
        Creates a SubRequest based on "self" and calls the execute method of the given fixturedef object. This will
        force the FixtureDef object to throw away any previous results and compute a new fixture value, which
        will be stored into the FixtureDef object itself.

        :param FixtureDef fixturedef:
        """
        # prepare a subrequest object before calling fixture function
        # (latter managed by fixturedef)
        argname = fixturedef.argname
        funcitem = self._pyfuncitem
        scope = fixturedef.scope
        try:
            param = funcitem.callspec.getparam(argname)
        except (AttributeError, ValueError):
            param = NOTSET
            param_index = 0
            if fixturedef.params is not None:
                frame = inspect.stack()[3]
                frameinfo = inspect.getframeinfo(frame[0])
                source_path = frameinfo.filename
                source_lineno = frameinfo.lineno
                source_path = py.path.local(source_path)
                if source_path.relto(funcitem.config.rootdir):
                    source_path = source_path.relto(funcitem.config.rootdir)
                msg = (
                    "The requested fixture has no parameter defined for the "
                    "current test.\n\nRequested fixture '{0}' defined in:\n{1}"
                    "\n\nRequested here:\n{2}:{3}".format(
                        fixturedef.argname,
                        getlocation(fixturedef.func, funcitem.config.rootdir),
                        source_path,
                        source_lineno,
                    )
                )
                fail(msg)
        else:
            # indices might not be set if old-style metafunc.addcall() was used
            param_index = funcitem.callspec.indices.get(argname, 0)
            # if a parametrize invocation set a scope it will override
            # the static scope defined with the fixture function
            paramscopenum = funcitem.callspec._arg2scopenum.get(argname)
            if paramscopenum is not None:
                scope = scopes[paramscopenum]

        subrequest = SubRequest(self, scope, param, param_index, fixturedef)

        # check if a higher-level scoped fixture accesses a lower level one
        subrequest._check_scope(argname, self.scope, scope)

        # clear sys.exc_info before invoking the fixture (python bug?)
        # if its not explicitly cleared it will leak into the call
        exc_clear()
        try:
            # call the fixture function
            fixturedef.execute(request=subrequest)
        finally:
            # if fixture function failed it might have registered finalizers
            self.session._setupstate.addfinalizer(functools.partial(fixturedef.finish, request=subrequest),
                                                  subrequest.node)

    def _check_scope(self, argname, invoking_scope, requested_scope):
        if argname == "request":
            return
        if scopemismatch(invoking_scope, requested_scope):
            # try to report something helpful
            lines = self._factorytraceback()
            fail("ScopeMismatch: You tried to access the %r scoped "
                 "fixture %r with a %r scoped request object, "
                 "involved factories\n%s" % (
                     (requested_scope, argname, invoking_scope, "\n".join(lines))),
                 pytrace=False)

    def _factorytraceback(self):
        lines = []
        for fixturedef in self._get_fixturestack():
            factory = fixturedef.func
            fs, lineno = getfslineno(factory)
            p = self._pyfuncitem.session.fspath.bestrelpath(fs)
            args = _format_args(factory)
            lines.append("%s:%d:  def %s%s" % (
                p, lineno, factory.__name__, args))
        return lines

    def _getscopeitem(self, scope):
        if scope == "function":
            # this might also be a non-function Item despite its attribute name
            return self._pyfuncitem
        node = get_scope_node(self._pyfuncitem, scope)
        if node is None and scope == "class":
            # fallback to function item itself
            node = self._pyfuncitem
        assert node, 'Could not obtain a node for scope "{}" for function {!r}'.format(scope, self._pyfuncitem)
        return node

    def __repr__(self):
        return "<FixtureRequest for %r>" % (self.node)


class SubRequest(FixtureRequest):
    """ a sub request for handling getting a fixture from a
    test function/fixture. """

    def __init__(self, request, scope, param, param_index, fixturedef):
        self._parent_request = request
        self.fixturename = fixturedef.argname
        if param is not NOTSET:
            self.param = param
        self.param_index = param_index
        self.scope = scope
        self._fixturedef = fixturedef
        self._pyfuncitem = request._pyfuncitem
        self._fixture_defs = request._fixture_defs
        self._arg2fixturedefs = request._arg2fixturedefs
        self._arg2index = request._arg2index
        self._fixturemanager = request._fixturemanager

    def __repr__(self):
        return "<SubRequest %r for %r>" % (self.fixturename, self._pyfuncitem)

    def addfinalizer(self, finalizer):
        self._fixturedef.addfinalizer(finalizer)


class ScopeMismatchError(Exception):
    """ A fixture function tries to use a different fixture function which
    which has a lower scope (e.g. a Session one calls a function one)
    """


scopes = "session module class function".split()
scopenum_function = scopes.index("function")


def scopemismatch(currentscope, newscope):
    return scopes.index(newscope) > scopes.index(currentscope)


def scope2index(scope, descr, where=None):
    """Look up the index of ``scope`` and raise a descriptive value error
    if not defined.
    """
    try:
        return scopes.index(scope)
    except ValueError:
        raise ValueError(
            "{0} {1}has an unsupported scope value '{2}'".format(
                descr, 'from {0} '.format(where) if where else '',
                scope)
        )


class FixtureLookupError(LookupError):
    """ could not return a requested Fixture (missing or invalid). """

    def __init__(self, argname, request, msg=None):
        self.argname = argname
        self.request = request
        self.fixturestack = request._get_fixturestack()
        self.msg = msg

    def formatrepr(self):
        tblines = []
        addline = tblines.append
        stack = [self.request._pyfuncitem.obj]
        stack.extend(map(lambda x: x.func, self.fixturestack))
        msg = self.msg
        if msg is not None:
            # the last fixture raise an error, let's present
            # it at the requesting side
            stack = stack[:-1]
        for function in stack:
            fspath, lineno = getfslineno(function)
            try:
                lines, _ = inspect.getsourcelines(get_real_func(function))
            except (IOError, IndexError, TypeError):
                error_msg = "file %s, line %s: source code not available"
                addline(error_msg % (fspath, lineno + 1))
            else:
                addline("file %s, line %s" % (fspath, lineno + 1))
                for i, line in enumerate(lines):
                    line = line.rstrip()
                    addline("  " + line)
                    if line.lstrip().startswith('def'):
                        break

        if msg is None:
            fm = self.request._fixturemanager
            available = []
            parentid = self.request._pyfuncitem.parent.nodeid
            for name, fixturedefs in fm._arg2fixturedefs.items():
                faclist = list(fm._matchfactories(fixturedefs, parentid))
                if faclist and name not in available:
                    available.append(name)
            msg = "fixture %r not found" % (self.argname,)
            msg += "\n available fixtures: %s" % (", ".join(sorted(available)),)
            msg += "\n use 'pytest --fixtures [testpath]' for help on them."

        return FixtureLookupErrorRepr(fspath, lineno, tblines, msg, self.argname)


class FixtureLookupErrorRepr(TerminalRepr):
    def __init__(self, filename, firstlineno, tblines, errorstring, argname):
        self.tblines = tblines
        self.errorstring = errorstring
        self.filename = filename
        self.firstlineno = firstlineno
        self.argname = argname

    def toterminal(self, tw):
        # tw.line("FixtureLookupError: %s" %(self.argname), red=True)
        for tbline in self.tblines:
            tw.line(tbline.rstrip())
        lines = self.errorstring.split("\n")
        if lines:
            tw.line('{0}       {1}'.format(FormattedExcinfo.fail_marker,
                                           lines[0].strip()), red=True)
            for line in lines[1:]:
                tw.line('{0}       {1}'.format(FormattedExcinfo.flow_marker,
                                               line.strip()), red=True)
        tw.line()
        tw.line("%s:%d" % (self.filename, self.firstlineno + 1))


def fail_fixturefunc(fixturefunc, msg):
    fs, lineno = getfslineno(fixturefunc)
    location = "%s:%s" % (fs, lineno + 1)
    source = _pytest._code.Source(fixturefunc)
    fail(msg + ":\n\n" + str(source.indent()) + "\n" + location,
         pytrace=False)


def call_fixture_func(fixturefunc, request, kwargs):
    yieldctx = is_generator(fixturefunc)
    if yieldctx:
        it = fixturefunc(**kwargs)
        res = next(it)

        def teardown():
            try:
                next(it)
            except StopIteration:
                pass
            else:
                fail_fixturefunc(fixturefunc,
                                 "yield_fixture function has more than one 'yield'")

        request.addfinalizer(teardown)
    else:
        res = fixturefunc(**kwargs)
    return res


class FixtureDef(object):
    """ A container for a factory definition. """

    def __init__(self, fixturemanager, baseid, argname, func, scope, params,
                 unittest=False, ids=None):
        self._fixturemanager = fixturemanager
        self.baseid = baseid or ''
        self.has_location = baseid is not None
        self.func = func
        self.argname = argname
        self.scope = scope
        self.scopenum = scope2index(
            scope or "function",
            descr='fixture {0}'.format(func.__name__),
            where=baseid
        )
        self.params = params
        self.argnames = getfuncargnames(func, is_method=unittest)
        self.unittest = unittest
        self.ids = ids
        self._finalizers = []

    def addfinalizer(self, finalizer):
        self._finalizers.append(finalizer)

    def finish(self, request):
        exceptions = []
        try:
            while self._finalizers:
                try:
                    func = self._finalizers.pop()
                    func()
                except:  # noqa
                    exceptions.append(sys.exc_info())
            if exceptions:
                e = exceptions[0]
                del exceptions  # ensure we don't keep all frames alive because of the traceback
                py.builtin._reraise(*e)

        finally:
            hook = self._fixturemanager.session.gethookproxy(request.node.fspath)
            hook.pytest_fixture_post_finalizer(fixturedef=self, request=request)
            # even if finalization fails, we invalidate
            # the cached fixture value and remove
            # all finalizers because they may be bound methods which will
            # keep instances alive
            if hasattr(self, "cached_result"):
                del self.cached_result
            self._finalizers = []

    def execute(self, request):
        # get required arguments and register our own finish()
        # with their finalization
        for argname in self.argnames:
            fixturedef = request._get_active_fixturedef(argname)
            if argname != "request":
                fixturedef.addfinalizer(functools.partial(self.finish, request=request))

        my_cache_key = request.param_index
        cached_result = getattr(self, "cached_result", None)
        if cached_result is not None:
            result, cache_key, err = cached_result
            if my_cache_key == cache_key:
                if err is not None:
                    py.builtin._reraise(*err)
                else:
                    return result
            # we have a previous but differently parametrized fixture instance
            # so we need to tear it down before creating a new one
            self.finish(request)
            assert not hasattr(self, "cached_result")

        hook = self._fixturemanager.session.gethookproxy(request.node.fspath)
        return hook.pytest_fixture_setup(fixturedef=self, request=request)

    def __repr__(self):
        return ("<FixtureDef name=%r scope=%r baseid=%r >" %
                (self.argname, self.scope, self.baseid))


def pytest_fixture_setup(fixturedef, request):
    """ Execution of fixture setup. """
    kwargs = {}
    for argname in fixturedef.argnames:
        fixdef = request._get_active_fixturedef(argname)
        result, arg_cache_key, exc = fixdef.cached_result
        request._check_scope(argname, request.scope, fixdef.scope)
        kwargs[argname] = result

    fixturefunc = fixturedef.func
    if fixturedef.unittest:
        if request.instance is not None:
            # bind the unbound method to the TestCase instance
            fixturefunc = fixturedef.func.__get__(request.instance)
    else:
        # the fixture function needs to be bound to the actual
        # request.instance so that code working with "fixturedef" behaves
        # as expected.
        if request.instance is not None:
            fixturefunc = getimfunc(fixturedef.func)
            if fixturefunc != fixturedef.func:
                fixturefunc = fixturefunc.__get__(request.instance)
    my_cache_key = request.param_index
    try:
        result = call_fixture_func(fixturefunc, request, kwargs)
    except TEST_OUTCOME:
        fixturedef.cached_result = (None, my_cache_key, sys.exc_info())
        raise
    fixturedef.cached_result = (result, my_cache_key, None)
    return result


def _ensure_immutable_ids(ids):
    if ids is None:
        return
    if callable(ids):
        return ids
    return tuple(ids)


@attr.s(frozen=True)
class FixtureFunctionMarker(object):
    scope = attr.ib()
    params = attr.ib(converter=attr.converters.optional(tuple))
    autouse = attr.ib(default=False)
    ids = attr.ib(default=None, converter=_ensure_immutable_ids)
    name = attr.ib(default=None)

    def __call__(self, function):
        if isclass(function):
            raise ValueError(
                "class fixtures not supported (may be in the future)")
        function._pytestfixturefunction = self
        return function


def fixture(scope="function", params=None, autouse=False, ids=None, name=None):
    """Decorator to mark a fixture factory function.

    This decorator can be used (with or without parameters) to define a
    fixture function.  The name of the fixture function can later be
    referenced to cause its invocation ahead of running tests: test
    modules or classes can use the pytest.mark.usefixtures(fixturename)
    marker.  Test functions can directly use fixture names as input
    arguments in which case the fixture instance returned from the fixture
    function will be injected.

    :arg scope: the scope for which this fixture is shared, one of
                "function" (default), "class", "module" or "session".

    :arg params: an optional list of parameters which will cause multiple
                invocations of the fixture function and all of the tests
                using it.

    :arg autouse: if True, the fixture func is activated for all tests that
                can see it.  If False (the default) then an explicit
                reference is needed to activate the fixture.

    :arg ids: list of string ids each corresponding to the params
                so that they are part of the test id. If no ids are provided
                they will be generated automatically from the params.

    :arg name: the name of the fixture. This defaults to the name of the
                decorated function. If a fixture is used in the same module in
                which it is defined, the function name of the fixture will be
                shadowed by the function arg that requests the fixture; one way
                to resolve this is to name the decorated function
                ``fixture_<fixturename>`` and then use
                ``@pytest.fixture(name='<fixturename>')``.

    Fixtures can optionally provide their values to test functions using a ``yield`` statement,
    instead of ``return``. In this case, the code block after the ``yield`` statement is executed
    as teardown code regardless of the test outcome. A fixture function must yield exactly once.
    """
    if callable(scope) and params is None and autouse is False:
        # direct decoration
        return FixtureFunctionMarker(
            "function", params, autouse, name=name)(scope)
    if params is not None and not isinstance(params, (list, tuple)):
        params = list(params)
    return FixtureFunctionMarker(scope, params, autouse, ids=ids, name=name)


def yield_fixture(scope="function", params=None, autouse=False, ids=None, name=None):
    """ (return a) decorator to mark a yield-fixture factory function.

    .. deprecated:: 3.0
        Use :py:func:`pytest.fixture` directly instead.
    """
    if callable(scope) and params is None and not autouse:
        # direct decoration
        return FixtureFunctionMarker(
            "function", params, autouse, ids=ids, name=name)(scope)
    else:
        return FixtureFunctionMarker(scope, params, autouse, ids=ids, name=name)


defaultfuncargprefixmarker = fixture()


@fixture(scope="session")
def pytestconfig(request):
    """Session-scoped fixture that returns the :class:`_pytest.config.Config` object.

    Example::

        def test_foo(pytestconfig):
            if pytestconfig.getoption("verbose"):
                ...

    """
    return request.config


class FixtureManager(object):
    """
    pytest fixtures definitions and information is stored and managed
    from this class.

    During collection fm.parsefactories() is called multiple times to parse
    fixture function definitions into FixtureDef objects and internal
    data structures.

    During collection of test functions, metafunc-mechanics instantiate
    a FuncFixtureInfo object which is cached per node/func-name.
    This FuncFixtureInfo object is later retrieved by Function nodes
    which themselves offer a fixturenames attribute.

    The FuncFixtureInfo object holds information about fixtures and FixtureDefs
    relevant for a particular function.  An initial list of fixtures is
    assembled like this:

    - ini-defined usefixtures
    - autouse-marked fixtures along the collection chain up from the function
    - usefixtures markers at module/class/function level
    - test function funcargs

    Subsequently the funcfixtureinfo.fixturenames attribute is computed
    as the closure of the fixtures needed to setup the initial fixtures,
    i. e. fixtures needed by fixture functions themselves are appended
    to the fixturenames list.

    Upon the test-setup phases all fixturenames are instantiated, retrieved
    by a lookup of their FuncFixtureInfo.
    """

    _argprefix = "pytest_funcarg__"
    FixtureLookupError = FixtureLookupError
    FixtureLookupErrorRepr = FixtureLookupErrorRepr

    def __init__(self, session):
        self.session = session
        self.config = session.config
        self._arg2fixturedefs = {}
        self._holderobjseen = set()
        self._arg2finish = {}
        self._nodeid_and_autousenames = [("", self.config.getini("usefixtures"))]
        session.config.pluginmanager.register(self, "funcmanage")

    def getfixtureinfo(self, node, func, cls, funcargs=True):
        if funcargs and not hasattr(node, "nofuncargs"):
            argnames = getfuncargnames(func, cls=cls)
        else:
            argnames = ()
        usefixtures = getattr(func, "usefixtures", None)
        initialnames = argnames
        if usefixtures is not None:
            initialnames = usefixtures.args + initialnames
        fm = node.session._fixturemanager
        names_closure, arg2fixturedefs = fm.getfixtureclosure(initialnames,
                                                              node)
        return FuncFixtureInfo(argnames, names_closure, arg2fixturedefs)

    def pytest_plugin_registered(self, plugin):
        nodeid = None
        try:
            p = py.path.local(plugin.__file__)
        except AttributeError:
            pass
        else:
            # construct the base nodeid which is later used to check
            # what fixtures are visible for particular tests (as denoted
            # by their test id)
            if p.basename.startswith("conftest.py"):
                nodeid = p.dirpath().relto(self.config.rootdir)
                if p.sep != nodes.SEP:
                    nodeid = nodeid.replace(p.sep, nodes.SEP)
        self.parsefactories(plugin, nodeid)

    def _getautousenames(self, nodeid):
        """ return a tuple of fixture names to be used. """
        autousenames = []
        for baseid, basenames in self._nodeid_and_autousenames:
            if nodeid.startswith(baseid):
                if baseid:
                    i = len(baseid)
                    nextchar = nodeid[i:i + 1]
                    if nextchar and nextchar not in ":/":
                        continue
                autousenames.extend(basenames)
        return autousenames

    def getfixtureclosure(self, fixturenames, parentnode):
        # collect the closure of all fixtures , starting with the given
        # fixturenames as the initial set.  As we have to visit all
        # factory definitions anyway, we also return a arg2fixturedefs
        # mapping so that the caller can reuse it and does not have
        # to re-discover fixturedefs again for each fixturename
        # (discovering matching fixtures for a given name/node is expensive)

        parentid = parentnode.nodeid
        fixturenames_closure = self._getautousenames(parentid)

        def merge(otherlist):
            for arg in otherlist:
                if arg not in fixturenames_closure:
                    fixturenames_closure.append(arg)

        merge(fixturenames)
        arg2fixturedefs = {}
        lastlen = -1
        while lastlen != len(fixturenames_closure):
            lastlen = len(fixturenames_closure)
            for argname in fixturenames_closure:
                if argname in arg2fixturedefs:
                    continue
                fixturedefs = self.getfixturedefs(argname, parentid)
                if fixturedefs:
                    arg2fixturedefs[argname] = fixturedefs
                    merge(fixturedefs[-1].argnames)

        def sort_by_scope(arg_name):
            try:
                fixturedefs = arg2fixturedefs[arg_name]
            except KeyError:
                return scopes.index('function')
            else:
                return fixturedefs[-1].scopenum

        fixturenames_closure.sort(key=sort_by_scope)
        return fixturenames_closure, arg2fixturedefs

    def pytest_generate_tests(self, metafunc):
        for argname in metafunc.fixturenames:
            faclist = metafunc._arg2fixturedefs.get(argname)
            if faclist:
                fixturedef = faclist[-1]
                if fixturedef.params is not None:
                    parametrize_func = getattr(metafunc.function, 'parametrize', None)
                    func_params = getattr(parametrize_func, 'args', [[None]])
                    func_kwargs = getattr(parametrize_func, 'kwargs', {})
                    # skip directly parametrized arguments
                    if "argnames" in func_kwargs:
                        argnames = parametrize_func.kwargs["argnames"]
                    else:
                        argnames = func_params[0]
                    if not isinstance(argnames, (tuple, list)):
                        argnames = [x.strip() for x in argnames.split(",") if x.strip()]
                    if argname not in func_params and argname not in argnames:
                        metafunc.parametrize(argname, fixturedef.params,
                                             indirect=True, scope=fixturedef.scope,
                                             ids=fixturedef.ids)
            else:
                continue  # will raise FixtureLookupError at setup time

    def pytest_collection_modifyitems(self, items):
        # separate parametrized setups
        items[:] = reorder_items(items)

    def parsefactories(self, node_or_obj, nodeid=NOTSET, unittest=False):
        if nodeid is not NOTSET:
            holderobj = node_or_obj
        else:
            holderobj = node_or_obj.obj
            nodeid = node_or_obj.nodeid
        if holderobj in self._holderobjseen:
            return
        self._holderobjseen.add(holderobj)
        autousenames = []
        for name in dir(holderobj):
            # The attribute can be an arbitrary descriptor, so the attribute
            # access below can raise. safe_getatt() ignores such exceptions.
            obj = safe_getattr(holderobj, name, None)
            # fixture functions have a pytest_funcarg__ prefix (pre-2.3 style)
            # or are "@pytest.fixture" marked
            marker = getfixturemarker(obj)
            if marker is None:
                if not name.startswith(self._argprefix):
                    continue
                if not callable(obj):
                    continue
                marker = defaultfuncargprefixmarker
                from _pytest import deprecated
                self.config.warn('C1', deprecated.FUNCARG_PREFIX.format(name=name), nodeid=nodeid)
                name = name[len(self._argprefix):]
            elif not isinstance(marker, FixtureFunctionMarker):
                # magic globals  with __getattr__ might have got us a wrong
                # fixture attribute
                continue
            else:
                if marker.name:
                    name = marker.name
                msg = 'fixtures cannot have "pytest_funcarg__" prefix ' \
                      'and be decorated with @pytest.fixture:\n%s' % name
                assert not name.startswith(self._argprefix), msg

            fixture_def = FixtureDef(self, nodeid, name, obj,
                                     marker.scope, marker.params,
                                     unittest=unittest, ids=marker.ids)

            faclist = self._arg2fixturedefs.setdefault(name, [])
            if fixture_def.has_location:
                faclist.append(fixture_def)
            else:
                # fixturedefs with no location are at the front
                # so this inserts the current fixturedef after the
                # existing fixturedefs from external plugins but
                # before the fixturedefs provided in conftests.
                i = len([f for f in faclist if not f.has_location])
                faclist.insert(i, fixture_def)
            if marker.autouse:
                autousenames.append(name)

        if autousenames:
            self._nodeid_and_autousenames.append((nodeid or '', autousenames))

    def getfixturedefs(self, argname, nodeid):
        """
        Gets a list of fixtures which are applicable to the given node id.

        :param str argname: name of the fixture to search for
        :param str nodeid: full node id of the requesting test.
        :return: list[FixtureDef]
        """
        try:
            fixturedefs = self._arg2fixturedefs[argname]
        except KeyError:
            return None
        else:
            return tuple(self._matchfactories(fixturedefs, nodeid))

    def _matchfactories(self, fixturedefs, nodeid):
        for fixturedef in fixturedefs:
            if nodes.ischildnode(fixturedef.baseid, nodeid):
                yield fixturedef
