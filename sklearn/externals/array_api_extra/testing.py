"""
Public testing utilities.

See also _lib._testing for additional private testing utilities.
"""

# https://github.com/scikit-learn/scikit-learn/pull/27910#issuecomment-2568023972
from __future__ import annotations

import contextlib
from collections.abc import Callable, Iterable, Iterator, Sequence
from functools import wraps
from types import ModuleType
from typing import TYPE_CHECKING, Any, TypeVar, cast

from ._lib._utils._compat import is_dask_namespace, is_jax_namespace

__all__ = ["lazy_xp_function", "patch_lazy_xp_functions"]

if TYPE_CHECKING:  # pragma: no cover
    # TODO move ParamSpec outside TYPE_CHECKING
    # depends on scikit-learn abandoning Python 3.9
    # https://github.com/scikit-learn/scikit-learn/pull/27910#issuecomment-2568023972
    from typing import ParamSpec

    import pytest
    from dask.typing import Graph, Key, SchedulerGetCallable
    from typing_extensions import override

    P = ParamSpec("P")
else:
    SchedulerGetCallable = object

    # Sphinx hacks
    class P:  # pylint: disable=missing-class-docstring
        args: tuple
        kwargs: dict

    def override(func: Callable[P, T]) -> Callable[P, T]:
        return func


T = TypeVar("T")

_ufuncs_tags: dict[object, dict[str, Any]] = {}  # type: ignore[explicit-any]


def lazy_xp_function(  # type: ignore[explicit-any]
    func: Callable[..., Any],
    *,
    allow_dask_compute: int = 0,
    jax_jit: bool = True,
    static_argnums: int | Sequence[int] | None = None,
    static_argnames: str | Iterable[str] | None = None,
) -> None:  # numpydoc ignore=GL07
    """
    Tag a function to be tested on lazy backends.

    Tag a function so that when any tests are executed with ``xp=jax.numpy`` the
    function is replaced with a jitted version of itself, and when it is executed with
    ``xp=dask.array`` the function will raise if it attempts to materialize the graph.
    This will be later expanded to provide test coverage for other lazy backends.

    In order for the tag to be effective, the test or a fixture must call
    :func:`patch_lazy_xp_functions`.

    Parameters
    ----------
    func : callable
        Function to be tested.
    allow_dask_compute : int, optional
        Number of times `func` is allowed to internally materialize the Dask graph. This
        is typically triggered by ``bool()``, ``float()``, or ``np.asarray()``.

        Set to 1 if you are aware that `func` converts the input parameters to numpy and
        want to let it do so at least for the time being, knowing that it is going to be
        extremely detrimental for performance.

        If a test needs values higher than 1 to pass, it is a canary that the conversion
        to numpy/bool/float is happening multiple times, which translates to multiple
        computations of the whole graph. Short of making the function fully lazy, you
        should at least add explicit calls to ``np.asarray()`` early in the function.
        *Note:* the counter of `allow_dask_compute` resets after each call to `func`, so
        a test function that invokes `func` multiple times should still work with this
        parameter set to 1.

        Default: 0, meaning that `func` must be fully lazy and never materialize the
        graph.
    jax_jit : bool, optional
        Set to True to replace `func` with ``jax.jit(func)`` after calling the
        :func:`patch_lazy_xp_functions` test helper with ``xp=jax.numpy``. Set to False
        if `func` is only compatible with eager (non-jitted) JAX. Default: True.
    static_argnums : int | Sequence[int], optional
        Passed to jax.jit. Positional arguments to treat as static (compile-time
        constant). Default: infer from `static_argnames` using
        `inspect.signature(func)`.
    static_argnames : str | Iterable[str], optional
        Passed to jax.jit. Named arguments to treat as static (compile-time constant).
        Default: infer from `static_argnums` using `inspect.signature(func)`.

    See Also
    --------
    patch_lazy_xp_functions : Companion function to call from the test or fixture.
    jax.jit : JAX function to compile a function for performance.

    Examples
    --------
    In ``test_mymodule.py``::

      from array_api_extra.testing import lazy_xp_function from mymodule import myfunc

      lazy_xp_function(myfunc)

      def test_myfunc(xp):
          a = xp.asarray([1, 2])
          # When xp=jax.numpy, this is the same as `b = jax.jit(myfunc)(a)`
          # When xp=dask.array, crash on compute() or persist()
          b = myfunc(a)

    Notes
    -----
    In order for this tag to be effective, the test function must be imported into the
    test module globals without its namespace; alternatively its namespace must be
    declared in a ``lazy_xp_modules`` list in the test module globals.

    Example 1::

      from mymodule import myfunc

      lazy_xp_function(myfunc)

      def test_myfunc(xp):
          x = myfunc(xp.asarray([1, 2]))

    Example 2::

      import mymodule

      lazy_xp_modules = [mymodule]
      lazy_xp_function(mymodule.myfunc)

      def test_myfunc(xp):
          x = mymodule.myfunc(xp.asarray([1, 2]))

    A test function can circumvent this monkey-patching system by using a namespace
    outside of the two above patterns. You need to sanitize your code to make sure this
    only happens intentionally.

    Example 1::

      import mymodule
      from mymodule import myfunc

      lazy_xp_function(myfunc)

      def test_myfunc(xp):
          a = xp.asarray([1, 2])
          b = myfunc(a)  # This is wrapped when xp=jax.numpy or xp=dask.array
          c = mymodule.myfunc(a)  # This is not

    Example 2::

      import mymodule

      class naked:
          myfunc = mymodule.myfunc

      lazy_xp_modules = [mymodule]
      lazy_xp_function(mymodule.myfunc)

      def test_myfunc(xp):
          a = xp.asarray([1, 2])
          b = mymodule.myfunc(a)  # This is wrapped when xp=jax.numpy or xp=dask.array
          c = naked.myfunc(a)  # This is not
    """
    tags = {
        "allow_dask_compute": allow_dask_compute,
        "jax_jit": jax_jit,
        "static_argnums": static_argnums,
        "static_argnames": static_argnames,
    }
    try:
        func._lazy_xp_function = tags  # type: ignore[attr-defined]  # pylint: disable=protected-access  # pyright: ignore[reportFunctionMemberAccess]
    except AttributeError:  # @cython.vectorize
        _ufuncs_tags[func] = tags


def patch_lazy_xp_functions(
    request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch, *, xp: ModuleType
) -> None:
    """
    Test lazy execution of functions tagged with :func:`lazy_xp_function`.

    If ``xp==jax.numpy``, search for all functions which have been tagged with
    :func:`lazy_xp_function` in the globals of the module that defines the current test,
    as well as in the ``lazy_xp_modules`` list in the globals of the same module,
    and wrap them with :func:`jax.jit`. Unwrap them at the end of the test.

    If ``xp==dask.array``, wrap the functions with a decorator that disables
    ``compute()`` and ``persist()`` and ensures that exceptions and warnings are raised
    eagerly.

    This function should be typically called by your library's `xp` fixture that runs
    tests on multiple backends::

        @pytest.fixture(params=[numpy, array_api_strict, jax.numpy, dask.array])
        def xp(request, monkeypatch):
            patch_lazy_xp_functions(request, monkeypatch, xp=request.param)
            return request.param

    but it can be otherwise be called by the test itself too.

    Parameters
    ----------
    request : pytest.FixtureRequest
        Pytest fixture, as acquired by the test itself or by one of its fixtures.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture, as acquired by the test itself or by one of its fixtures.
    xp : array_namespace
        Array namespace to be tested.

    See Also
    --------
    lazy_xp_function : Tag a function to be tested on lazy backends.
    pytest.FixtureRequest : `request` test function parameter.
    """
    mod = cast(ModuleType, request.module)
    mods = [mod, *cast(list[ModuleType], getattr(mod, "lazy_xp_modules", []))]

    def iter_tagged() -> (  # type: ignore[explicit-any]
        Iterator[tuple[ModuleType, str, Callable[..., Any], dict[str, Any]]]
    ):
        for mod in mods:
            for name, func in mod.__dict__.items():
                tags: dict[str, Any] | None = None  # type: ignore[explicit-any]
                with contextlib.suppress(AttributeError):
                    tags = func._lazy_xp_function  # pylint: disable=protected-access
                if tags is None:
                    with contextlib.suppress(KeyError, TypeError):
                        tags = _ufuncs_tags[func]
                if tags is not None:
                    yield mod, name, func, tags

    if is_dask_namespace(xp):
        for mod, name, func, tags in iter_tagged():
            n = tags["allow_dask_compute"]
            wrapped = _dask_wrap(func, n)
            monkeypatch.setattr(mod, name, wrapped)

    elif is_jax_namespace(xp):
        import jax

        for mod, name, func, tags in iter_tagged():
            if tags["jax_jit"]:
                # suppress unused-ignore to run mypy in -e lint as well as -e dev
                wrapped = cast(  # type: ignore[explicit-any]
                    Callable[..., Any],
                    jax.jit(
                        func,
                        static_argnums=tags["static_argnums"],
                        static_argnames=tags["static_argnames"],
                    ),
                )
                monkeypatch.setattr(mod, name, wrapped)


class CountingDaskScheduler(SchedulerGetCallable):
    """
    Dask scheduler that counts how many times `dask.compute` is called.

    If the number of times exceeds 'max_count', it raises an error.
    This is a wrapper around Dask's own 'synchronous' scheduler.

    Parameters
    ----------
    max_count : int
        Maximum number of allowed calls to `dask.compute`.
    msg : str
        Assertion to raise when the count exceeds `max_count`.
    """

    count: int
    max_count: int
    msg: str

    def __init__(self, max_count: int, msg: str):  # numpydoc ignore=GL08
        self.count = 0
        self.max_count = max_count
        self.msg = msg

    @override
    def __call__(self, dsk: Graph, keys: Sequence[Key] | Key, **kwargs: Any) -> Any:  # type: ignore[decorated-any,explicit-any] # numpydoc ignore=GL08
        import dask

        self.count += 1
        # This should yield a nice traceback to the
        # offending line in the user's code
        assert self.count <= self.max_count, self.msg

        return dask.get(dsk, keys, **kwargs)  # type: ignore[attr-defined,no-untyped-call] # pyright: ignore[reportPrivateImportUsage]


def _dask_wrap(
    func: Callable[P, T], n: int
) -> Callable[P, T]:  # numpydoc ignore=PR01,RT01
    """
    Wrap `func` to raise if it attempts to call `dask.compute` more than `n` times.

    After the function returns, materialize the graph in order to re-raise exceptions.
    """
    import dask

    func_name = getattr(func, "__name__", str(func))
    n_str = f"only up to {n}" if n else "no"
    msg = (
        f"Called `dask.compute()` or `dask.persist()` {n + 1} times, "
        f"but {n_str} calls are allowed. Set "
        f"`lazy_xp_function({func_name}, allow_dask_compute={n + 1})` "
        "to allow for more (but note that this will harm performance). "
    )

    @wraps(func)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:  # numpydoc ignore=GL08
        scheduler = CountingDaskScheduler(n, msg)
        with dask.config.set({"scheduler": scheduler}):  # pyright: ignore[reportPrivateImportUsage]
            out = func(*args, **kwargs)

        # Block until the graph materializes and reraise exceptions. This allows
        # `pytest.raises` and `pytest.warns` to work as expected. Note that this would
        # not work on scheduler='distributed', as it would not block.
        return dask.persist(out, scheduler="threads")[0]  # type: ignore[attr-defined,no-untyped-call,func-returns-value,index]  # pyright: ignore[reportPrivateImportUsage]

    return wrapper
