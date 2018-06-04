import sys
import types
import toolz
from toolz.compatibility import import_module


class TlzLoader(object):
    """ Finds and loads ``tlz`` modules when added to sys.meta_path"""
    def __init__(self):
        self.always_from_toolz = set([
            toolz.pipe,
        ])

    def _load_toolz(self, fullname):
        rv = {}
        package, dot, submodules = fullname.partition('.')
        try:
            module_name = ''.join(['cytoolz', dot, submodules])
            rv['cytoolz'] = import_module(module_name)
        except ImportError:
            pass
        try:
            module_name = ''.join(['toolz', dot, submodules])
            rv['toolz'] = import_module(module_name)
        except ImportError:
            pass
        if not rv:
            raise ImportError(fullname)
        return rv

    def find_module(self, fullname, path=None):  # pragma: py3 no cover
        package, dot, submodules = fullname.partition('.')
        if package == 'tlz':
            return self

    def load_module(self, fullname):  # pragma: py3 no cover
        if fullname in sys.modules:  # pragma: no cover
            return sys.modules[fullname]
        spec = TlzSpec(fullname, self)
        module = self.create_module(spec)
        sys.modules[fullname] = module
        self.exec_module(module)
        return module

    def find_spec(self, fullname, path, target=None):  # pragma: no cover
        package, dot, submodules = fullname.partition('.')
        if package == 'tlz':
            return TlzSpec(fullname, self)

    def create_module(self, spec):
        return types.ModuleType(spec.name)

    def exec_module(self, module):
        toolz_mods = self._load_toolz(module.__name__)
        fast_mod = toolz_mods.get('cytoolz') or toolz_mods['toolz']
        slow_mod = toolz_mods.get('toolz') or toolz_mods['cytoolz']
        module.__dict__.update(toolz.merge(fast_mod.__dict__, module.__dict__))
        package = fast_mod.__package__
        if package is not None:
            package, dot, submodules = package.partition('.')
            module.__package__ = ''.join(['tlz', dot, submodules])
        if not module.__doc__:
            module.__doc__ = fast_mod.__doc__

        # show file from toolz during introspection
        module.__file__ = slow_mod.__file__

        for k, v in fast_mod.__dict__.items():
            tv = slow_mod.__dict__.get(k)
            try:
                hash(tv)
            except TypeError:
                tv = None
            if tv in self.always_from_toolz:
                module.__dict__[k] = tv
            elif (
                isinstance(v, types.ModuleType)
                and v.__package__ == fast_mod.__name__
            ):
                package, dot, submodules = v.__name__.partition('.')
                module_name = ''.join(['tlz', dot, submodules])
                submodule = import_module(module_name)
                module.__dict__[k] = submodule


class TlzSpec(object):
    def __init__(self, name, loader):
        self.name = name
        self.loader = loader
        self.origin = None
        self.submodule_search_locations = []
        self.loader_state = None
        self.cached = None
        self.parent = None
        self.has_location = False


tlz_loader = TlzLoader()
sys.meta_path.append(tlz_loader)
tlz_loader.exec_module(sys.modules['tlz'])
