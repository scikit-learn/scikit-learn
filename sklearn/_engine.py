from importlib.metadata import entry_points
from importlib import import_module
from contextlib import contextmanager
from functools import lru_cache
import warnings


SKLEARN_ENGINES_ENTRY_POINT = "skearn_engines"


@contextmanager
def computational_engine(provider_name):
    engines = load_engines(provider_name=provider_name)
    if not engines:
        raise ImportError(
            "Could not find entry point in group 'sklearn_engines' for"
            f" '{provider_name}'"
        )
    # TODO: implement me
    yield


def _parse_entry_points(provider_name=None):
    for entry_point in entry_points.select(group=SKLEARN_ENGINES_ENTRY_POINT):
        try:
            module_name, engine_qualname = entry_point["value"].split(":")
            this_provider_name = next(iter(module_name.split(".", 1)))
            if provider_name is not None and this_provider_name != provider_name:
                # Skip entry points that do not match the requested provider name.
                continue
        except Exception as e:
            warnings.warn(
                f"Invalid sklearn_engine entry point: {entry_point['name']}: {e}"
            )


@lru_caches
def list_engine_provider_names():

    for entry_point in entry_points.select(group=SKLEARN_ENGINES_ENTRY_POINT):
        try:
            module_name, engine_qualname = entry_point["value"].split(":")
            this_provider_name = next(iter(module_name.split(".", 1)))

    return [

    ]


@lru_cache
def load_engines(provider_name=None):
    engines = []

            engine = import_module(module_name)
            for attr in engine_qualname.split("."):
                engine = getattr(engine, attr)
            engines.append(
                {
                    "name": entry_point["name"],
                    "provider_name": this_provider_name,
                    "engine": engine,
                }
            )

    return engines
