"""Find all objects reachable from a root object."""

from collections.abc import Iterable
import weakref
import types

from typing import List, Dict, Iterator, Tuple, Mapping
from typing_extensions import Final

method_descriptor_type: Final = type(object.__dir__)
method_wrapper_type: Final = type(object().__ne__)
wrapper_descriptor_type: Final = type(object.__ne__)

FUNCTION_TYPES: Final = (
    types.BuiltinFunctionType,
    types.FunctionType,
    types.MethodType,
    method_descriptor_type,
    wrapper_descriptor_type,
    method_wrapper_type,
)

ATTR_BLACKLIST: Final = {
    '__doc__',
    '__name__',
    '__class__',
    '__dict__',
}

# Instances of these types can't have references to other objects
ATOMIC_TYPE_BLACKLIST: Final = {
    bool,
    int,
    float,
    str,
    type(None),
    object,
}

# Don't look at most attributes of these types
COLLECTION_TYPE_BLACKLIST: Final = {
    list,
    set,
    dict,
    tuple,
}

# Don't return these objects
TYPE_BLACKLIST: Final = {
    weakref.ReferenceType,
}


def isproperty(o: object, attr: str) -> bool:
    return isinstance(getattr(type(o), attr, None), property)


def get_edge_candidates(o: object) -> Iterator[Tuple[object, object]]:
    # use getattr because mypyc expects dict, not mappingproxy
    if '__getattribute__' in getattr(type(o), '__dict__'):  # noqa
        return
    if type(o) not in COLLECTION_TYPE_BLACKLIST:
        for attr in dir(o):
            try:
                if attr not in ATTR_BLACKLIST and hasattr(o, attr) and not isproperty(o, attr):
                    e = getattr(o, attr)
                    if not type(e) in ATOMIC_TYPE_BLACKLIST:
                        yield attr, e
            except AssertionError:
                pass
    if isinstance(o, Mapping):
        for k, v in o.items():
            yield k, v
    elif isinstance(o, Iterable) and not isinstance(o, str):
        for i, e in enumerate(o):
            yield i, e


def get_edges(o: object) -> Iterator[Tuple[object, object]]:
    for s, e in get_edge_candidates(o):
        if (isinstance(e, FUNCTION_TYPES)):
            # We don't want to collect methods, but do want to collect values
            # in closures and self pointers to other objects

            if hasattr(e, '__closure__'):
                yield (s, '__closure__'), e.__closure__  # type: ignore
            if hasattr(e, '__self__'):
                se = e.__self__  # type: ignore
                if se is not o and se is not type(o) and hasattr(s, '__self__'):
                    yield s.__self__, se  # type: ignore
        else:
            if not type(e) in TYPE_BLACKLIST:
                yield s, e


def get_reachable_graph(root: object) -> Tuple[Dict[int, object],
                                               Dict[int, Tuple[int, object]]]:
    parents = {}
    seen = {id(root): root}
    worklist = [root]
    while worklist:
        o = worklist.pop()
        for s, e in get_edges(o):
            if id(e) in seen:
                continue
            parents[id(e)] = (id(o), s)
            seen[id(e)] = e
            worklist.append(e)

    return seen, parents


def get_path(o: object,
             seen: Dict[int, object],
             parents: Dict[int, Tuple[int, object]]) -> List[Tuple[object, object]]:
    path = []
    while id(o) in parents:
        pid, attr = parents[id(o)]
        o = seen[pid]
        path.append((attr, o))
    path.reverse()
    return path
