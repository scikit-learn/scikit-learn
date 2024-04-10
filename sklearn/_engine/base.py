import inspect
import warnings
from functools import lru_cache, wraps
from importlib import import_module
from importlib.metadata import entry_points

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
        engine_entry_points = all_entry_points.get(SKLEARN_ENGINES_ENTRY_POINT, ())
    for entry_point in engine_entry_points:
        try:
            spec = _parse_entry_point(entry_point)
            if provider_names is not None and spec.provider_name not in provider_names:
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


def _get_engine_classes(engine_name, provider_names, engine_specs, default):
    specs_by_provider = {}
    for spec in engine_specs:
        if spec.name != engine_name:
            continue
        specs_by_provider.setdefault(spec.provider_name, spec)

    for provider_name in provider_names:
        if inspect.isclass(provider_name):
            # The provider name is actually a ready-to-go engine class.
            # Instead of a made up string to name this ad-hoc provider
            # we use the class itself. This mirrors what the user used
            # when they set the config (ad-hoc class or string naming
            # a provider).
            engine_class = provider_name
            if getattr(engine_class, "engine_name", None) != engine_name:
                continue
            yield engine_class, engine_class

        spec = specs_by_provider.get(provider_name)
        if spec is not None:
            yield spec.provider_name, spec.get_engine_class()

    yield "default", default


def get_engine_classes(engine_name, default, verbose=False):
    """Find all possible providers of `engine_name`.

    Provider candidates are found based on parsing entrypoint definitions that
    match the name of enabled engine providers, as well as, ad-hoc providers
    in the form of engine classes in the list of enabled engine providers.

    Parameters
    ----------
    engine_name : str
        The name of the algorithm for which to find engine classes.

    default : class
        The default engine class to use if no other provider is found.

    verbose : bool, default=False
        If True, print the name of the engine classes that are tried.

    Yields
    ------
    provider : str or class
        The "name" of each matching provider. The "name" corresponds to the
        entry in the `engine_provider` configuration. It can be a string or a
        class for programmatically registered ad-hoc providers.

    engine_class :
        The engine class that implements the algorithm for the given provider.
    """
    provider_names = get_config()["engine_provider"]

    if not provider_names:
        yield "default", default
        return

    engine_specs = _parse_entry_points(
        provider_names=tuple(
            [name for name in provider_names if not inspect.isclass(name)]
        )
    )
    for provider, engine_class in _get_engine_classes(
        engine_name=engine_name,
        provider_names=provider_names,
        engine_specs=engine_specs,
        default=default,
    ):
        if verbose:
            print(
                f"trying engine {engine_class.__module__}.{engine_class.__qualname__}."
            )
        yield provider, engine_class


def convert_attributes(method):
    """Convert estimator attributes after calling the decorated method.

    The attributes of an estimator can be stored in "engine native" types
    (default) or "scikit-learn native" types. This decorator will call the
    engine's conversion function when needed. Use this decorator on methods
    that set estimator attributes.
    """

    @wraps(method)
    def wrapper(self, *args, **kwargs):
        r = method(self, *args, **kwargs)
        convert_attributes = get_config()["engine_attributes"]

        if convert_attributes == "sklearn_types":
            engine = self._engine_class
            for name, value in vars(self).items():
                # All attributes are passed to the engine, which can
                # either convert the value (engine specific types) or
                # return it as is (native Python types)
                converted = engine.convert_to_sklearn_types(name, value)
                setattr(self, name, converted)

            # No matter which engine was used to fit, after the attribute
            # conversion to the sklearn native types the default engine
            # is used.
            self._engine_class = self._default_engine
            self._engine_provider = "default"

        return r

    return wrapper
