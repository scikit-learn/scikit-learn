import inspect
import warnings
from .callers import _multicall, HookCallError, _Result, _legacymulticall

__version__ = '0.6.0'

__all__ = ["PluginManager", "PluginValidationError", "HookCallError",
           "HookspecMarker", "HookimplMarker"]


class PluginValidationError(Exception):
    """ plugin failed validation. """


class HookspecMarker(object):
    """ Decorator helper class for marking functions as hook specifications.

    You can instantiate it with a project_name to get a decorator.
    Calling PluginManager.add_hookspecs later will discover all marked functions
    if the PluginManager uses the same project_name.
    """

    def __init__(self, project_name):
        self.project_name = project_name

    def __call__(self, function=None, firstresult=False, historic=False):
        """ if passed a function, directly sets attributes on the function
        which will make it discoverable to add_hookspecs().  If passed no
        function, returns a decorator which can be applied to a function
        later using the attributes supplied.

        If firstresult is True the 1:N hook call (N being the number of registered
        hook implementation functions) will stop at I<=N when the I'th function
        returns a non-None result.

        If historic is True calls to a hook will be memorized and replayed
        on later registered plugins.

        """
        def setattr_hookspec_opts(func):
            if historic and firstresult:
                raise ValueError("cannot have a historic firstresult hook")
            setattr(func, self.project_name + "_spec",
                    dict(firstresult=firstresult, historic=historic))
            return func

        if function is not None:
            return setattr_hookspec_opts(function)
        else:
            return setattr_hookspec_opts


class HookimplMarker(object):
    """ Decorator helper class for marking functions as hook implementations.

    You can instantiate with a project_name to get a decorator.
    Calling PluginManager.register later will discover all marked functions
    if the PluginManager uses the same project_name.
    """
    def __init__(self, project_name):
        self.project_name = project_name

    def __call__(self, function=None, hookwrapper=False, optionalhook=False,
                 tryfirst=False, trylast=False):

        """ if passed a function, directly sets attributes on the function
        which will make it discoverable to register().  If passed no function,
        returns a decorator which can be applied to a function later using
        the attributes supplied.

        If optionalhook is True a missing matching hook specification will not result
        in an error (by default it is an error if no matching spec is found).

        If tryfirst is True this hook implementation will run as early as possible
        in the chain of N hook implementations for a specfication.

        If trylast is True this hook implementation will run as late as possible
        in the chain of N hook implementations.

        If hookwrapper is True the hook implementations needs to execute exactly
        one "yield".  The code before the yield is run early before any non-hookwrapper
        function is run.  The code after the yield is run after all non-hookwrapper
        function have run.  The yield receives a ``_Result`` object representing
        the exception or result outcome of the inner calls (including other hookwrapper
        calls).

        """
        def setattr_hookimpl_opts(func):
            setattr(func, self.project_name + "_impl",
                    dict(hookwrapper=hookwrapper, optionalhook=optionalhook,
                         tryfirst=tryfirst, trylast=trylast))
            return func

        if function is None:
            return setattr_hookimpl_opts
        else:
            return setattr_hookimpl_opts(function)


def normalize_hookimpl_opts(opts):
    opts.setdefault("tryfirst", False)
    opts.setdefault("trylast", False)
    opts.setdefault("hookwrapper", False)
    opts.setdefault("optionalhook", False)


class _TagTracer(object):
    def __init__(self):
        self._tag2proc = {}
        self.writer = None
        self.indent = 0

    def get(self, name):
        return _TagTracerSub(self, (name,))

    def format_message(self, tags, args):
        if isinstance(args[-1], dict):
            extra = args[-1]
            args = args[:-1]
        else:
            extra = {}

        content = " ".join(map(str, args))
        indent = "  " * self.indent

        lines = [
            "%s%s [%s]\n" % (indent, content, ":".join(tags))
        ]

        for name, value in extra.items():
            lines.append("%s    %s: %s\n" % (indent, name, value))
        return lines

    def processmessage(self, tags, args):
        if self.writer is not None and args:
            lines = self.format_message(tags, args)
            self.writer(''.join(lines))
        try:
            self._tag2proc[tags](tags, args)
        except KeyError:
            pass

    def setwriter(self, writer):
        self.writer = writer

    def setprocessor(self, tags, processor):
        if isinstance(tags, str):
            tags = tuple(tags.split(":"))
        else:
            assert isinstance(tags, tuple)
        self._tag2proc[tags] = processor


class _TagTracerSub(object):
    def __init__(self, root, tags):
        self.root = root
        self.tags = tags

    def __call__(self, *args):
        self.root.processmessage(self.tags, args)

    def setmyprocessor(self, processor):
        self.root.setprocessor(self.tags, processor)

    def get(self, name):
        return self.__class__(self.root, self.tags + (name,))


class _TracedHookExecution(object):
    def __init__(self, pluginmanager, before, after):
        self.pluginmanager = pluginmanager
        self.before = before
        self.after = after
        self.oldcall = pluginmanager._inner_hookexec
        assert not isinstance(self.oldcall, _TracedHookExecution)
        self.pluginmanager._inner_hookexec = self

    def __call__(self, hook, hook_impls, kwargs):
        self.before(hook.name, hook_impls, kwargs)
        outcome = _Result.from_call(lambda: self.oldcall(hook, hook_impls, kwargs))
        self.after(outcome, hook.name, hook_impls, kwargs)
        return outcome.get_result()

    def undo(self):
        self.pluginmanager._inner_hookexec = self.oldcall


class PluginManager(object):
    """ Core Pluginmanager class which manages registration
    of plugin objects and 1:N hook calling.

    You can register new hooks by calling ``add_hookspec(module_or_class)``.
    You can register plugin objects (which contain hooks) by calling
    ``register(plugin)``.  The Pluginmanager is initialized with a
    prefix that is searched for in the names of the dict of registered
    plugin objects.  An optional excludefunc allows to blacklist names which
    are not considered as hooks despite a matching prefix.

    For debugging purposes you can call ``enable_tracing()``
    which will subsequently send debug information to the trace helper.
    """

    def __init__(self, project_name, implprefix=None):
        """ if implprefix is given implementation functions
        will be recognized if their name matches the implprefix. """
        self.project_name = project_name
        self._name2plugin = {}
        self._plugin2hookcallers = {}
        self._plugin_distinfo = []
        self.trace = _TagTracer().get("pluginmanage")
        self.hook = _HookRelay(self.trace.root.get("hook"))
        self._implprefix = implprefix
        self._inner_hookexec = lambda hook, methods, kwargs: \
            hook.multicall(
                methods, kwargs,
                firstresult=hook.spec_opts.get('firstresult'),
            )

    def _hookexec(self, hook, methods, kwargs):
        # called from all hookcaller instances.
        # enable_tracing will set its own wrapping function at self._inner_hookexec
        return self._inner_hookexec(hook, methods, kwargs)

    def register(self, plugin, name=None):
        """ Register a plugin and return its canonical name or None if the name
        is blocked from registering.  Raise a ValueError if the plugin is already
        registered. """
        plugin_name = name or self.get_canonical_name(plugin)

        if plugin_name in self._name2plugin or plugin in self._plugin2hookcallers:
            if self._name2plugin.get(plugin_name, -1) is None:
                return  # blocked plugin, return None to indicate no registration
            raise ValueError("Plugin already registered: %s=%s\n%s" %
                             (plugin_name, plugin, self._name2plugin))

        # XXX if an error happens we should make sure no state has been
        # changed at point of return
        self._name2plugin[plugin_name] = plugin

        # register matching hook implementations of the plugin
        self._plugin2hookcallers[plugin] = hookcallers = []
        for name in dir(plugin):
            hookimpl_opts = self.parse_hookimpl_opts(plugin, name)
            if hookimpl_opts is not None:
                normalize_hookimpl_opts(hookimpl_opts)
                method = getattr(plugin, name)
                hookimpl = HookImpl(plugin, plugin_name, method, hookimpl_opts)
                hook = getattr(self.hook, name, None)
                if hook is None:
                    hook = _HookCaller(name, self._hookexec)
                    setattr(self.hook, name, hook)
                elif hook.has_spec():
                    self._verify_hook(hook, hookimpl)
                    hook._maybe_apply_history(hookimpl)
                hook._add_hookimpl(hookimpl)
                hookcallers.append(hook)
        return plugin_name

    def parse_hookimpl_opts(self, plugin, name):
        method = getattr(plugin, name)
        if not inspect.isroutine(method):
            return
        try:
            res = getattr(method, self.project_name + "_impl", None)
        except Exception:
            res = {}
        if res is not None and not isinstance(res, dict):
            # false positive
            res = None
        elif res is None and self._implprefix and name.startswith(self._implprefix):
            res = {}
        return res

    def unregister(self, plugin=None, name=None):
        """ unregister a plugin object and all its contained hook implementations
        from internal data structures. """
        if name is None:
            assert plugin is not None, "one of name or plugin needs to be specified"
            name = self.get_name(plugin)

        if plugin is None:
            plugin = self.get_plugin(name)

        # if self._name2plugin[name] == None registration was blocked: ignore
        if self._name2plugin.get(name):
            del self._name2plugin[name]

        for hookcaller in self._plugin2hookcallers.pop(plugin, []):
            hookcaller._remove_plugin(plugin)

        return plugin

    def set_blocked(self, name):
        """ block registrations of the given name, unregister if already registered. """
        self.unregister(name=name)
        self._name2plugin[name] = None

    def is_blocked(self, name):
        """ return True if the name blogs registering plugins of that name. """
        return name in self._name2plugin and self._name2plugin[name] is None

    def add_hookspecs(self, module_or_class):
        """ add new hook specifications defined in the given module_or_class.
        Functions are recognized if they have been decorated accordingly. """
        names = []
        for name in dir(module_or_class):
            spec_opts = self.parse_hookspec_opts(module_or_class, name)
            if spec_opts is not None:
                hc = getattr(self.hook, name, None)
                if hc is None:
                    hc = _HookCaller(name, self._hookexec, module_or_class, spec_opts)
                    setattr(self.hook, name, hc)
                else:
                    # plugins registered this hook without knowing the spec
                    hc.set_specification(module_or_class, spec_opts)
                    for hookfunction in (hc._wrappers + hc._nonwrappers):
                        self._verify_hook(hc, hookfunction)
                names.append(name)

        if not names:
            raise ValueError("did not find any %r hooks in %r" %
                             (self.project_name, module_or_class))

    def parse_hookspec_opts(self, module_or_class, name):
        method = getattr(module_or_class, name)
        return getattr(method, self.project_name + "_spec", None)

    def get_plugins(self):
        """ return the set of registered plugins. """
        return set(self._plugin2hookcallers)

    def is_registered(self, plugin):
        """ Return True if the plugin is already registered. """
        return plugin in self._plugin2hookcallers

    def get_canonical_name(self, plugin):
        """ Return canonical name for a plugin object. Note that a plugin
        may be registered under a different name which was specified
        by the caller of register(plugin, name). To obtain the name
        of an registered plugin use ``get_name(plugin)`` instead."""
        return getattr(plugin, "__name__", None) or str(id(plugin))

    def get_plugin(self, name):
        """ Return a plugin or None for the given name. """
        return self._name2plugin.get(name)

    def has_plugin(self, name):
        """ Return True if a plugin with the given name is registered. """
        return self.get_plugin(name) is not None

    def get_name(self, plugin):
        """ Return name for registered plugin or None if not registered. """
        for name, val in self._name2plugin.items():
            if plugin == val:
                return name

    def _verify_hook(self, hook, hookimpl):
        if hook.is_historic() and hookimpl.hookwrapper:
            raise PluginValidationError(
                "Plugin %r\nhook %r\nhistoric incompatible to hookwrapper" %
                (hookimpl.plugin_name, hook.name))

        # positional arg checking
        notinspec = set(hookimpl.argnames) - set(hook.argnames)
        if notinspec:
            raise PluginValidationError(
                "Plugin %r for hook %r\nhookimpl definition: %s\n"
                "Argument(s) %s are declared in the hookimpl but "
                "can not be found in the hookspec" %
                (hookimpl.plugin_name, hook.name,
                 _formatdef(hookimpl.function), notinspec)
            )

    def check_pending(self):
        """ Verify that all hooks which have not been verified against
        a hook specification are optional, otherwise raise PluginValidationError"""
        for name in self.hook.__dict__:
            if name[0] != "_":
                hook = getattr(self.hook, name)
                if not hook.has_spec():
                    for hookimpl in (hook._wrappers + hook._nonwrappers):
                        if not hookimpl.optionalhook:
                            raise PluginValidationError(
                                "unknown hook %r in plugin %r" %
                                (name, hookimpl.plugin))

    def load_setuptools_entrypoints(self, entrypoint_name):
        """ Load modules from querying the specified setuptools entrypoint name.
        Return the number of loaded plugins. """
        from pkg_resources import (iter_entry_points, DistributionNotFound,
                                   VersionConflict)
        for ep in iter_entry_points(entrypoint_name):
            # is the plugin registered or blocked?
            if self.get_plugin(ep.name) or self.is_blocked(ep.name):
                continue
            try:
                plugin = ep.load()
            except DistributionNotFound:
                continue
            except VersionConflict as e:
                raise PluginValidationError(
                    "Plugin %r could not be loaded: %s!" % (ep.name, e))
            self.register(plugin, name=ep.name)
            self._plugin_distinfo.append((plugin, ep.dist))
        return len(self._plugin_distinfo)

    def list_plugin_distinfo(self):
        """ return list of distinfo/plugin tuples for all setuptools registered
        plugins. """
        return list(self._plugin_distinfo)

    def list_name_plugin(self):
        """ return list of name/plugin pairs. """
        return list(self._name2plugin.items())

    def get_hookcallers(self, plugin):
        """ get all hook callers for the specified plugin. """
        return self._plugin2hookcallers.get(plugin)

    def add_hookcall_monitoring(self, before, after):
        """ add before/after tracing functions for all hooks
        and return an undo function which, when called,
        will remove the added tracers.

        ``before(hook_name, hook_impls, kwargs)`` will be called ahead
        of all hook calls and receive a hookcaller instance, a list
        of HookImpl instances and the keyword arguments for the hook call.

        ``after(outcome, hook_name, hook_impls, kwargs)`` receives the
        same arguments as ``before`` but also a :py:class:`_Result`` object
        which represents the result of the overall hook call.
        """
        return _TracedHookExecution(self, before, after).undo

    def enable_tracing(self):
        """ enable tracing of hook calls and return an undo function. """
        hooktrace = self.hook._trace

        def before(hook_name, methods, kwargs):
            hooktrace.root.indent += 1
            hooktrace(hook_name, kwargs)

        def after(outcome, hook_name, methods, kwargs):
            if outcome.excinfo is None:
                hooktrace("finish", hook_name, "-->", outcome.get_result())
            hooktrace.root.indent -= 1

        return self.add_hookcall_monitoring(before, after)

    def subset_hook_caller(self, name, remove_plugins):
        """ Return a new _HookCaller instance for the named method
        which manages calls to all registered plugins except the
        ones from remove_plugins. """
        orig = getattr(self.hook, name)
        plugins_to_remove = [plug for plug in remove_plugins if hasattr(plug, name)]
        if plugins_to_remove:
            hc = _HookCaller(orig.name, orig._hookexec, orig._specmodule_or_class,
                             orig.spec_opts)
            for hookimpl in (orig._wrappers + orig._nonwrappers):
                plugin = hookimpl.plugin
                if plugin not in plugins_to_remove:
                    hc._add_hookimpl(hookimpl)
                    # we also keep track of this hook caller so it
                    # gets properly removed on plugin unregistration
                    self._plugin2hookcallers.setdefault(plugin, []).append(hc)
            return hc
        return orig


def varnames(func):
    """Return tuple of positional and keywrord argument names for a function,
    method, class or callable.

    In case of a class, its ``__init__`` method is considered.
    For methods the ``self`` parameter is not included.
    """
    cache = getattr(func, "__dict__", {})
    try:
        return cache["_varnames"]
    except KeyError:
        pass

    if inspect.isclass(func):
        try:
            func = func.__init__
        except AttributeError:
            return (), ()
    elif not inspect.isroutine(func):  # callable object?
        try:
            func = getattr(func, '__call__', func)
        except Exception:
            return ()

    try:  # func MUST be a function or method here or we won't parse any args
        spec = _getargspec(func)
    except TypeError:
        return (), ()

    args, defaults = tuple(spec.args), spec.defaults
    if defaults:
        index = -len(defaults)
        args, defaults = args[:index], tuple(args[index:])
    else:
        defaults = ()

    # strip any implicit instance arg
    if args:
        if inspect.ismethod(func) or (
            '.' in getattr(func, '__qualname__', ()) and args[0] == 'self'
        ):
            args = args[1:]

    assert "self" not in args  # best naming practises check?
    try:
        cache["_varnames"] = args, defaults
    except TypeError:
        pass
    return args, defaults


class _HookRelay(object):
    """ hook holder object for performing 1:N hook calls where N is the number
    of registered plugins.

    """

    def __init__(self, trace):
        self._trace = trace


class _HookCaller(object):
    def __init__(self, name, hook_execute, specmodule_or_class=None,
                 spec_opts=None):
        self.name = name
        self._wrappers = []
        self._nonwrappers = []
        self._hookexec = hook_execute
        self._specmodule_or_class = None
        self.argnames = None
        self.kwargnames = None
        self.multicall = _multicall
        self.spec_opts = spec_opts or {}
        if specmodule_or_class is not None:
            self.set_specification(specmodule_or_class, spec_opts)

    def has_spec(self):
        return self._specmodule_or_class is not None

    def set_specification(self, specmodule_or_class, spec_opts):
        assert not self.has_spec()
        self._specmodule_or_class = specmodule_or_class
        specfunc = getattr(specmodule_or_class, self.name)
        # get spec arg signature
        argnames, self.kwargnames = varnames(specfunc)
        self.argnames = ["__multicall__"] + list(argnames)
        self.spec_opts.update(spec_opts)
        if spec_opts.get("historic"):
            self._call_history = []

    def is_historic(self):
        return hasattr(self, "_call_history")

    def _remove_plugin(self, plugin):
        def remove(wrappers):
            for i, method in enumerate(wrappers):
                if method.plugin == plugin:
                    del wrappers[i]
                    return True
        if remove(self._wrappers) is None:
            if remove(self._nonwrappers) is None:
                raise ValueError("plugin %r not found" % (plugin,))

    def _add_hookimpl(self, hookimpl):
        """A an implementation to the callback chain.
        """
        if hookimpl.hookwrapper:
            methods = self._wrappers
        else:
            methods = self._nonwrappers

        if hookimpl.trylast:
            methods.insert(0, hookimpl)
        elif hookimpl.tryfirst:
            methods.append(hookimpl)
        else:
            # find last non-tryfirst method
            i = len(methods) - 1
            while i >= 0 and methods[i].tryfirst:
                i -= 1
            methods.insert(i + 1, hookimpl)

        if '__multicall__' in hookimpl.argnames:
            warnings.warn(
                "Support for __multicall__ is now deprecated and will be"
                "removed in an upcoming release.",
                DeprecationWarning
            )
            self.multicall = _legacymulticall

    def __repr__(self):
        return "<_HookCaller %r>" % (self.name,)

    def __call__(self, *args, **kwargs):
        if args:
            raise TypeError("hook calling supports only keyword arguments")
        assert not self.is_historic()
        if self.argnames:
            notincall = set(self.argnames) - set(['__multicall__']) - set(
                kwargs.keys())
            if notincall:
                warnings.warn(
                    "Argument(s) {} which are declared in the hookspec "
                    "can not be found in this hook call"
                    .format(tuple(notincall)),
                    stacklevel=2,
                )
        return self._hookexec(self, self._nonwrappers + self._wrappers, kwargs)

    def call_historic(self, proc=None, kwargs=None):
        """ call the hook with given ``kwargs`` for all registered plugins and
        for all plugins which will be registered afterwards.

        If ``proc`` is not None it will be called for for each non-None result
        obtained from a hook implementation.
        """
        self._call_history.append((kwargs or {}, proc))
        # historizing hooks don't return results
        res = self._hookexec(self, self._nonwrappers + self._wrappers, kwargs)
        for x in res or []:
            proc(x)

    def call_extra(self, methods, kwargs):
        """ Call the hook with some additional temporarily participating
        methods using the specified kwargs as call parameters. """
        old = list(self._nonwrappers), list(self._wrappers)
        for method in methods:
            opts = dict(hookwrapper=False, trylast=False, tryfirst=False)
            hookimpl = HookImpl(None, "<temp>", method, opts)
            self._add_hookimpl(hookimpl)
        try:
            return self(**kwargs)
        finally:
            self._nonwrappers, self._wrappers = old

    def _maybe_apply_history(self, method):
        """Apply call history to a new hookimpl if it is marked as historic.
        """
        if self.is_historic():
            for kwargs, proc in self._call_history:
                res = self._hookexec(self, [method], kwargs)
                if res and proc is not None:
                    proc(res[0])


class HookImpl(object):
    def __init__(self, plugin, plugin_name, function, hook_impl_opts):
        self.function = function
        self.argnames, self.kwargnames = varnames(self.function)
        self.plugin = plugin
        self.opts = hook_impl_opts
        self.plugin_name = plugin_name
        self.__dict__.update(hook_impl_opts)


if hasattr(inspect, 'getfullargspec'):
    def _getargspec(func):
        return inspect.getfullargspec(func)
else:
    def _getargspec(func):
        return inspect.getargspec(func)


if hasattr(inspect, 'signature'):
    def _formatdef(func):
        return "%s%s" % (
            func.__name__,
            str(inspect.signature(func))
        )
else:
    def _formatdef(func):
        return "%s%s" % (
            func.__name__,
            inspect.formatargspec(*inspect.getargspec(func))
        )
