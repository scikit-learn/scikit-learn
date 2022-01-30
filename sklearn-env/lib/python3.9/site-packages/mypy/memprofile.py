"""Utility for dumping memory usage stats.

This is tailored to mypy and knows (a little) about which list objects are
owned by particular AST nodes, etc.
"""

from collections import defaultdict
import gc
import sys
from typing import List, Dict, Iterable, Tuple, cast

from mypy.nodes import FakeInfo, Node
from mypy.types import Type
from mypy.util import get_class_descriptors


def collect_memory_stats() -> Tuple[Dict[str, int],
                                    Dict[str, int]]:
    """Return stats about memory use.

    Return a tuple with these items:
      - Dict from object kind to number of instances of that kind
      - Dict from object kind to total bytes used by all instances of that kind
    """
    objs = gc.get_objects()
    find_recursive_objects(objs)

    inferred = {}
    for obj in objs:
        if type(obj) is FakeInfo:
            # Processing these would cause a crash.
            continue
        n = type(obj).__name__
        if hasattr(obj, '__dict__'):
            # Keep track of which class a particular __dict__ is associated with.
            inferred[id(obj.__dict__)] = '%s (__dict__)' % n
        if isinstance(obj, (Node, Type)):  # type: ignore
            if hasattr(obj, '__dict__'):
                for x in obj.__dict__.values():
                    if isinstance(x, list):
                        # Keep track of which node a list is associated with.
                        inferred[id(x)] = '%s (list)' % n
                    if isinstance(x, tuple):
                        # Keep track of which node a list is associated with.
                        inferred[id(x)] = '%s (tuple)' % n

            for k in get_class_descriptors(type(obj)):
                x = getattr(obj, k, None)
                if isinstance(x, list):
                    inferred[id(x)] = '%s (list)' % n
                if isinstance(x, tuple):
                    inferred[id(x)] = '%s (tuple)' % n

    freqs: Dict[str, int] = {}
    memuse: Dict[str, int] = {}
    for obj in objs:
        if id(obj) in inferred:
            name = inferred[id(obj)]
        else:
            name = type(obj).__name__
        freqs[name] = freqs.get(name, 0) + 1
        memuse[name] = memuse.get(name, 0) + sys.getsizeof(obj)

    return freqs, memuse


def print_memory_profile(run_gc: bool = True) -> None:
    if not sys.platform.startswith('win'):
        import resource
        system_memuse = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    else:
        system_memuse = -1  # TODO: Support this on Windows
    if run_gc:
        gc.collect()
    freqs, memuse = collect_memory_stats()
    print('%7s  %7s  %7s  %s' % ('Freq', 'Size(k)', 'AvgSize', 'Type'))
    print('-------------------------------------------')
    totalmem = 0
    i = 0
    for n, mem in sorted(memuse.items(), key=lambda x: -x[1]):
        f = freqs[n]
        if i < 50:
            print('%7d  %7d  %7.0f  %s' % (f, mem // 1024, mem / f, n))
        i += 1
        totalmem += mem
    print()
    print('Mem usage RSS   ', system_memuse // 1024)
    print('Total reachable ', totalmem // 1024)


def find_recursive_objects(objs: List[object]) -> None:
    """Find additional objects referenced by objs and append them to objs.

    We use this since gc.get_objects() does not return objects without pointers
    in them such as strings.
    """
    seen = set(id(o) for o in objs)

    def visit(o: object) -> None:
        if id(o) not in seen:
            objs.append(o)
            seen.add(id(o))

    for obj in objs[:]:
        if type(obj) is FakeInfo:
            # Processing these would cause a crash.
            continue
        if type(obj) in (dict, defaultdict):
            for key, val in cast(Dict[object, object], obj).items():
                visit(key)
                visit(val)
        if type(obj) in (list, tuple, set):
            for x in cast(Iterable[object], obj):
                visit(x)
        if hasattr(obj, '__slots__'):
            for base in type.mro(type(obj)):
                for slot in getattr(base, '__slots__', ()):
                    if hasattr(obj, slot):
                        visit(getattr(obj, slot))
