from dataclasses import dataclass
from importlib.metadata import entry_points
from importlib import import_module
from contextlib import contextmanager
from functools import lru_cache
import warnings


SKLEARN_ENGINES_ENTRY_POINT = "skearn_engines"


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


@contextmanager
def computational_engine(provider_names):
    if isinstance(provider_names, str):
        provider_names = [provider_names]
    # TODO: implement me and complain if no entry point can be found with the given provider name
    parsed_entry_points = _parse_entry_points(provider_names=provider_names)
    if len(parsed_entry_points) == 0:
        raise RuntimeError()
    yield


def _parse_entry_point(entry_point):
    module_name, engine_qualname = entry_point["value"].split(":")
    provider_name = next(iter(module_name.split(".", 1)))
    return EngineSpec(entry_point["name"], provider_name, module_name, engine_qualname)


@lru_cache
def _parse_entry_points(provider_names=None):
    specs = []
    for entry_point in entry_points().select(group=SKLEARN_ENGINES_ENTRY_POINT):
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
                f"Invalid sklearn_engine entry point {entry_point['name']} "
                f"with value {entry_point['value']}: {e}"
            )
    return specs


def list_engine_provider_names():
    """Find the list of sklearn_engine provider names

    This function only inspects the metadata and should trigger any module import.
    """
    return sorted({spec.provider_name for spec in _parse_entry_points()})
