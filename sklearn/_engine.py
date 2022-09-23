from importlib.metadata import entry_points
from importlib import import_module
from functools import lru_cache
import warnings

from sklearn._config import get_config

SKLEARN_ENGINES_ENTRY_POINT = "sklearn_engines"


class EngineSpec:

    __slots__ = ["name", "provider_name", "module_name", "engine_qualname"]

    def __init__(self, name, provider_name, module_name, engine_qualname):
        self.name = name
        self.provider_name = provider_name
        self.module_name = module_name
        self.engine_qualname = engine_qualname

    def get_engine_class(self):
        engine = import_module(self.module_name)
        for attr in self.engine_qualname.split("."):
            engine = getattr(engine, attr)
        return engine


def _parse_entry_point(entry_point):
    module_name, engine_qualname = entry_point.value.split(":")
    provider_name = next(iter(module_name.split(".", 1)))
    return EngineSpec(entry_point.name, provider_name, module_name, engine_qualname)


@lru_cache
def _parse_entry_points(provider_names=None):
    specs = []
    all_entry_points = entry_points()
    if hasattr(all_entry_points, "select"):
        engine_entry_points = all_entry_points.select(group=SKLEARN_ENGINES_ENTRY_POINT)
    else:
        engine_entry_points = all_entry_points[SKLEARN_ENGINES_ENTRY_POINT]
    for entry_point in engine_entry_points:
        try:
            spec = _parse_entry_point(entry_point)
            if provider_names is not None and spec.provider_name in provider_names:
                # Skip entry points that do not match the requested provider names.
                continue
            specs.append(spec)
        except Exception as e:
            # Do not raise an exception in case an invalid package has been
            # installed in the same Python env as scikit-learn: just warn and
            # skip.
            warnings.warn(
                f"Invalid {SKLEARN_ENGINES_ENTRY_POINT} entry point"
                f" {entry_point.name} with value {entry_point.value}: {e}"
            )
    if provider_names is not None:
        observed_provider_names = {spec.provider_name for spec in specs}
        missing_providers = set(provider_names) - observed_provider_names
        if missing_providers:
            raise RuntimeError(
                "Could not find any provider for the"
                f" {SKLEARN_ENGINES_ENTRY_POINT} entry point with name(s):"
                f" {', '.join(repr(p) for p in sorted(missing_providers))}"
            )
    return specs


def list_engine_provider_names():
    """Find the list of sklearn_engine provider names

    This function only inspects the metadata and should trigger any module import.
    """
    return sorted({spec.provider_name for spec in _parse_entry_points()})


def _get_engine_class(engine_name, provider_names, engine_specs, default=None):
    specs_by_provider = {}
    for spec in engine_specs:
        if spec.name != engine_name:
            continue
        specs_by_provider.setdefault(spec.provider_name, spec)

    for provider_name in provider_names:
        spec = specs_by_provider.get(provider_name)
        if spec is not None:
            # XXX: should we return an instance or the class itself?
            return spec.get_engine_class()

    return default


def get_engine_class(engine_name, default=None):
    provider_names = get_config()["engine_provider"]
    if isinstance(provider_names, str):
        provider_names = (provider_names,)
    elif not isinstance(provider_names, tuple):
        # Make sure the provider names are a tuple to make it possible for the
        # lru cache to hash them.
        provider_names = tuple(provider_names)
    if not provider_names:
        return default
    engine_specs = _parse_entry_points(provider_names=provider_names)
    return _get_engine_class(
        engine_name=engine_name,
        provider_names=provider_names,
        engine_specs=engine_specs,
        default=default,
    )
