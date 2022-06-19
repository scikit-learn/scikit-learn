import itertools
import functools
import contextlib

from setuptools.extern.packaging.requirements import Requirement
from setuptools.extern.packaging.version import Version
from setuptools.extern.more_itertools import always_iterable
from setuptools.extern.jaraco.context import suppress
from setuptools.extern.jaraco.functools import apply

from ._compat import metadata, repair_extras


def resolve(req: Requirement) -> metadata.Distribution:
    """
    Resolve the requirement to its distribution.

    Ignore exception detail for Python 3.9 compatibility.

    >>> resolve(Requirement('pytest<3'))  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    importlib.metadata.PackageNotFoundError: No package metadata was found for pytest<3
    """
    dist = metadata.distribution(req.name)
    if not req.specifier.contains(Version(dist.version), prereleases=True):
        raise metadata.PackageNotFoundError(str(req))
    dist.extras = req.extras  # type: ignore
    return dist


@apply(bool)
@suppress(metadata.PackageNotFoundError)
def is_satisfied(req: Requirement):
    return resolve(req)


unsatisfied = functools.partial(itertools.filterfalse, is_satisfied)


class NullMarker:
    @classmethod
    def wrap(cls, req: Requirement):
        return req.marker or cls()

    def evaluate(self, *args, **kwargs):
        return True


def find_direct_dependencies(dist, extras=None):
    """
    Find direct, declared dependencies for dist.
    """
    simple = (
        req
        for req in map(Requirement, always_iterable(dist.requires))
        if NullMarker.wrap(req).evaluate(dict(extra=None))
    )
    extra_deps = (
        req
        for req in map(Requirement, always_iterable(dist.requires))
        for extra in always_iterable(getattr(dist, "extras", extras))
        if NullMarker.wrap(req).evaluate(dict(extra=extra))
    )
    return itertools.chain(simple, extra_deps)


def traverse(items, visit):
    """
    Given an iterable of items, traverse the items.

    For each item, visit is called to return any additional items
    to include in the traversal.
    """
    while True:
        try:
            item = next(items)
        except StopIteration:
            return
        yield item
        items = itertools.chain(items, visit(item))


def find_req_dependencies(req):
    with contextlib.suppress(metadata.PackageNotFoundError):
        dist = resolve(req)
        yield from find_direct_dependencies(dist)


def find_dependencies(dist, extras=None):
    """
    Find all reachable dependencies for dist.

    dist is an importlib.metadata.Distribution (or similar).
    TODO: create a suitable protocol for type hint.

    >>> deps = find_dependencies(resolve(Requirement('nspektr')))
    >>> all(isinstance(dep, Requirement) for dep in deps)
    True
    >>> not any('pytest' in str(dep) for dep in deps)
    True
    >>> test_deps = find_dependencies(resolve(Requirement('nspektr[testing]')))
    >>> any('pytest' in str(dep) for dep in test_deps)
    True
    """

    def visit(req, seen=set()):
        if req in seen:
            return ()
        seen.add(req)
        return find_req_dependencies(req)

    return traverse(find_direct_dependencies(dist, extras), visit)


class Unresolved(Exception):
    def __iter__(self):
        return iter(self.args[0])


def missing(ep):
    """
    Generate the unresolved dependencies (if any) of ep.
    """
    return unsatisfied(find_dependencies(ep.dist, repair_extras(ep.extras)))


def check(ep):
    """
    >>> ep, = metadata.entry_points(group='console_scripts', name='pip')
    >>> check(ep)
    >>> dist = metadata.distribution('nspektr')

    Since 'docs' extras are not installed, requesting them should fail.

    >>> ep = metadata.EntryPoint(
    ...     group=None, name=None, value='nspektr [docs]')._for(dist)
    >>> check(ep)
    Traceback (most recent call last):
    ...
    nspektr.Unresolved: [...]
    """
    missed = list(missing(ep))
    if missed:
        raise Unresolved(missed)
