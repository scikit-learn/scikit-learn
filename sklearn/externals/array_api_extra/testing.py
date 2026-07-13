"""
Public testing utilities.

See also _lib._testing for additional private testing utilities.
"""

from __future__ import annotations

import contextlib
import enum
import math
import warnings
from collections.abc import Callable, Generator, Iterator, Sequence
from functools import update_wrapper, wraps
from inspect import getattr_static
from types import FunctionType, ModuleType
from typing import TYPE_CHECKING, Any, ParamSpec, TypeVar, cast

from ._lib._utils._compat import (
    array_namespace,
    is_array_api_strict_namespace,
    is_cupy_namespace,
    is_dask_namespace,
    is_jax_namespace,
    is_numpy_namespace,
    is_pydata_sparse_namespace,
    is_torch_array,
    is_torch_namespace,
    to_device,
)
from ._lib._utils._helpers import jax_autojit, pickle_flatten, pickle_unflatten
from ._lib._utils._typing import Array, Device

__all__ = [
    "assert_close",
    "assert_close_nulp",
    "assert_equal",
    "assert_less",
    "lazy_xp_function",
    "patch_lazy_xp_functions",
]

if TYPE_CHECKING:  # pragma: no cover
    # TODO import override from typing (requires Python >=3.12)
    import numpy as np
    import pytest
    from dask.typing import Graph, Key, SchedulerGetCallable
    from typing_extensions import override

else:
    # Sphinx hacks
    SchedulerGetCallable = object

    def override(func):
        return func


P = ParamSpec("P")
T = TypeVar("T")

_ufuncs_tags: dict[object, dict[str, Any]] = {}


class Deprecated(enum.Enum):
    """Unique type for deprecated parameters."""

    DEPRECATED = 1


DEPRECATED = Deprecated.DEPRECATED


def _clone_function(  # numpydoc ignore=PR01,RT01
    f: Callable[..., Any],
) -> Callable[..., Any]:
    """Return a clone of an existing function."""
    f_new = FunctionType(
        f.__code__,
        f.__globals__,
        name=f.__name__,
        argdefs=f.__defaults__,
        closure=f.__closure__,
    )
    f_new.__kwdefaults__ = f.__kwdefaults__
    return update_wrapper(f_new, f)


def lazy_xp_function(
    func: Callable[..., Any] | tuple[type, str],
    *,
    allow_dask_compute: bool | int = False,
    jax_jit: bool = True,
    static_argnums: Deprecated = DEPRECATED,
    static_argnames: Deprecated = DEPRECATED,
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
    func : callable | tuple[type, str]
        Function to be tested, or a tuple containing an (uninstantiated) class and a
        method name to specify a class method to be tested.
    allow_dask_compute : bool | int, optional
        Whether `func` is allowed to internally materialize the Dask graph, or maximum
        number of times it is allowed to do so. This is typically triggered by
        ``bool()``, ``float()``, or ``np.asarray()``.

        Set to 1 if you are aware that `func` converts the input parameters to NumPy and
        want to let it do so at least for the time being, knowing that it is going to be
        extremely detrimental for performance.

        If a test needs values higher than 1 to pass, it is a canary that the conversion
        to NumPy/bool/float is happening multiple times, which translates to multiple
        computations of the whole graph. Short of making the function fully lazy, you
        should at least add explicit calls to ``np.asarray()`` early in the function.
        *Note:* the counter of `allow_dask_compute` resets after each call to `func`, so
        a test function that invokes `func` multiple times should still work with this
        parameter set to 1.

        Set to True to allow `func` to materialize the graph an unlimited number
        of times.

        Default: False, meaning that `func` must be fully lazy and never materialize the
        graph.
    jax_jit : bool, optional
        Set to True to replace `func` with a smart variant of ``jax.jit(func)`` after
        calling the :func:`patch_lazy_xp_functions` test helper with ``xp=jax.numpy``.
        This is the default behaviour.
        Set to False if `func` is only compatible with eager (non-jitted) JAX.

        Unlike with vanilla ``jax.jit``, all arguments and return types that are not JAX
        arrays are treated as static; the function can accept and return arbitrary
        wrappers around JAX arrays. This difference is because, in real life, most users
        won't wrap the function directly with ``jax.jit`` but rather they will use it
        within their own code, which is itself then wrapped by ``jax.jit``, and
        internally consume the function's outputs.

        In other words, the pattern that is being tested is::

            >>> @jax.jit
            ... def user_func(x):
            ...     y = user_prepares_inputs(x)
            ...     z = func(y, some_static_arg=True)
            ...     return user_consumes(z)

        Default: True.
    static_argnums : Deprecated
        Deprecated; ignored.
    static_argnames : Deprecated
        Deprecated; ignored.

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
          # When xp=jax.numpy, this is similar to `b = jax.jit(myfunc)(a)`
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
    if static_argnums is not DEPRECATED or static_argnames is not DEPRECATED:
        warnings.warn(
            (
                "The `static_argnums` and `static_argnames` parameters are deprecated "
                "and ignored. They will be removed in a future version."
            ),
            DeprecationWarning,
            stacklevel=2,
        )
    tags: dict[str, bool | int | type] = {
        "allow_dask_compute": allow_dask_compute,
        "jax_jit": jax_jit,
    }

    if isinstance(func, tuple):
        # Replace the method with a clone before adding tags
        # to avoid adding unwanted tags to a parent method when
        # the method was inherited from a parent class.
        # Note: can't just accept an unbound method `cls.method_name` because in
        # case of inheritance it would be impossible to attribute it to the child class.
        # This also makes it so tagged methods will appear in their class's ``__dict__``
        # and thus findable by ``iter_tagged_modules`` below.
        cls, method_name = func
        # The method might be a staticmethod or classmethod so we need to do a dance
        # to ensure that this is preserved.
        raw_attr = getattr_static(cls, method_name)
        method = getattr(cls, method_name)
        if isinstance(raw_attr, classmethod):
            method = method.__func__
        cloned_method = _clone_function(method)

        method_to_set: Any
        if isinstance(raw_attr, staticmethod):
            method_to_set = staticmethod(cloned_method)
        elif isinstance(raw_attr, classmethod):
            method_to_set = classmethod(cloned_method)
        else:
            method_to_set = cloned_method

        setattr(cls, method_name, method_to_set)
        f = getattr(cls, method_name)
        if isinstance(raw_attr, classmethod):
            f = f.__func__
        # Annotate that cls owns this method so we can check that later.
        tags["owner"] = cls
    else:
        f = func

    try:
        f._lazy_xp_function = tags  # pylint: disable=protected-access # pyright: ignore[reportFunctionMemberAccess] # pyrefly: ignore[missing-attribute]
    except AttributeError:  # @cython.vectorize
        _ufuncs_tags[f] = tags


def patch_lazy_xp_functions(
    request: pytest.FixtureRequest,
    monkeypatch: pytest.MonkeyPatch | None = None,
    *,
    xp: ModuleType,
) -> contextlib.AbstractContextManager[None]:
    """
    Test lazy execution of functions tagged with :func:`lazy_xp_function`.

    If ``xp==jax.numpy``, search for all functions and methods which have been tagged
    with :func:`lazy_xp_function` in the globals of the module that defines the current
    test, as well as in the ``lazy_xp_modules`` list in the globals of the same module,
    and wrap them with :func:`jax.jit`.
    Unwrap them at the end of the test.

    If ``xp==dask.array``, wrap the functions with a decorator that disables
    ``compute()`` and ``persist()`` and ensures that exceptions and warnings are raised
    eagerly.

    This function should be typically called by your library's `xp` fixture that runs
    tests on multiple backends::

        @pytest.fixture(params=[
            numpy,
            array_api_strict,
            pytest.param(jax.numpy, marks=pytest.mark.thread_unsafe),
            pytest.param(dask.array, marks=pytest.mark.thread_unsafe),
        ])
        def xp(request):
            with patch_lazy_xp_functions(request, xp=request.param):
                yield request.param

    but it can be otherwise be called by the test itself too.

    Parameters
    ----------
    request : pytest.FixtureRequest
        Pytest fixture, as acquired by the test itself or by one of its fixtures.
    monkeypatch : pytest.MonkeyPatch
        Deprecated.
    xp : array_namespace
        Array namespace to be tested.

    Returns
    -------
    contextlib.AbstractContextManager
        Testing context manager.

    See Also
    --------
    lazy_xp_function : Tag a function to be tested on lazy backends.
    pytest.FixtureRequest : `request` test function parameter.

    Notes
    -----
    This context manager monkey-patches modules and as such is thread unsafe
    on Dask and JAX. If you run your test suite with
    `pytest-run-parallel <https://github.com/Quansight-Labs/pytest-run-parallel/>`_,
    you should mark these backends with ``@pytest.mark.thread_unsafe``, as shown in
    the example above.
    """
    mod = cast(ModuleType, request.module)
    search_targets: list[ModuleType | type] = [
        mod,
        *cast(list[ModuleType], getattr(mod, "lazy_xp_modules", [])),
    ]
    # Also search for classes within the above modules which have had lazy_xp_function
    # applied to methods through ``lazy_xp_function((cls, method_name))`` syntax.
    # We might end up adding classes incidentally imported into modules, so using a
    # set here to cut down on potential redundancy.
    classes: set[type] = set()
    for target in search_targets:
        for obj in target.__dict__.values():
            if isinstance(obj, type):
                classes.add(obj)
    search_targets.extend(classes)

    to_revert: list[tuple[ModuleType | type, str, object]] = []

    def temp_setattr(  # numpydoc ignore=PR01
        target: ModuleType | type, name: str, func: object
    ) -> None:
        """
        Temporary setattr.

        Variant of monkeypatch.setattr, which allows monkey-patching only selected
        parameters of a test so that pytest-run-parallel can run on the remainder.
        """
        assert hasattr(target, name)
        # Need getattr_static because the attr could be a staticmethod or other
        # descriptor and we don't want that to be stripped away.
        original = getattr_static(target, name)
        to_revert.append((target, name, original))
        setattr(target, name, func)

    if monkeypatch is not None:
        warnings.warn(
            (
                "The `monkeypatch` parameter is deprecated and will be removed in a "
                "future version. "
                "Use `patch_lazy_xp_function` as a context manager instead."
            ),
            DeprecationWarning,
            stacklevel=2,
        )
        # Enable using patch_lazy_xp_function not as a context manager
        temp_setattr = monkeypatch.setattr  # type: ignore[assignment]  # pyright: ignore[reportAssignmentType]

    def iter_tagged() -> Iterator[
        tuple[ModuleType | type, str, Any, Callable[..., Any], dict[str, Any]]
    ]:  # numpydoc ignore=GL08
        for target in search_targets:
            for name, attr in target.__dict__.items():
                # attr might be a staticmethod or classmethod. If so we need
                # to peel it back and wrap the underlying function and later
                # make sure not to accidentally replace it with a regular
                # method.
                func: Any = (
                    attr.__func__
                    if isinstance(attr, (staticmethod, classmethod))
                    else attr
                )
                tags: dict[str, Any] | None = None
                with contextlib.suppress(AttributeError):
                    tags = func._lazy_xp_function  # pylint: disable=protected-access
                if tags is None:
                    with contextlib.suppress(KeyError, TypeError):
                        tags = _ufuncs_tags[func]
                if tags is not None:
                    if isinstance(target, type) and tags.get("owner") is not target:
                        # There's a common pattern to wrap functions in namespace
                        # classes to bypass lazy_xp_function like this:
                        #
                        # class naked:
                        #     myfunc = mymodule.myfunc
                        #
                        # To ensure this still works when checking for tags in
                        # attributes of classes, ensure that target is the actual
                        # owning class where func was defined.
                        continue
                    # put attr, and func in the outputs so we can later tell
                    # if this was a staticmethod or classmethod.
                    yield target, name, attr, func, tags

    wrapped: Any
    if is_dask_namespace(xp):
        for target, name, attr, func, tags in iter_tagged():
            n = tags["allow_dask_compute"]
            if n is True:
                n = 1_000_000
            elif n is False:
                n = 0
            wrapped = _dask_wrap(func, n)
            # If we're dealing with a staticmethod or classmethod, make
            # sure things stay that way.
            if isinstance(attr, staticmethod):
                wrapped = staticmethod(wrapped)
            elif isinstance(attr, classmethod):
                wrapped = classmethod(wrapped)
            temp_setattr(target, name, wrapped)

    elif is_jax_namespace(xp):
        for target, name, attr, func, tags in iter_tagged():
            if tags["jax_jit"]:
                wrapped = jax_autojit(func)
                # If we're dealing with a staticmethod or classmethod, make
                # sure things stay that way.
                if isinstance(attr, staticmethod):
                    wrapped = staticmethod(wrapped)
                elif isinstance(attr, classmethod):
                    wrapped = classmethod(wrapped)
                temp_setattr(target, name, wrapped)

    # We can't just decorate patch_lazy_xp_functions with
    # @contextlib.contextmanager because it would not work with the
    # deprecated monkeypatch when not used as a context manager.
    @contextlib.contextmanager
    def revert_on_exit() -> Generator[None]:  # numpydoc ignore=GL08
        try:
            yield
        finally:
            for target, name, orig_func in to_revert:
                setattr(target, name, orig_func)

    return revert_on_exit()


class _CountingDaskScheduler(SchedulerGetCallable):
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

    def __init__(self, max_count: int, msg: str) -> None:  # numpydoc ignore=GL08
        self.count = 0
        self.max_count = max_count
        self.msg = msg

    @override
    def __call__(
        self, dsk: Graph, keys: Sequence[Key] | Key, **kwargs: Any
    ) -> Any:  # numpydoc ignore=GL08
        import dask

        self.count += 1
        # This should yield a nice traceback to the
        # offending line in the user's code
        assert self.count <= self.max_count, self.msg

        return dask.get(dsk, keys, **kwargs)  # type: ignore[attr-defined]  # pyright: ignore[reportPrivateImportUsage]


def _dask_wrap(
    func: Callable[P, T], n: int
) -> Callable[P, T]:  # numpydoc ignore=PR01,RT01
    """
    Wrap `func` to raise if it attempts to call `dask.compute` more than `n` times.

    After the function returns, materialize the graph in order to re-raise exceptions.
    """
    import dask
    import dask.array as da

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
        scheduler = _CountingDaskScheduler(n, msg)
        with dask.config.set({"scheduler": scheduler}):  # pyright: ignore[reportPrivateImportUsage]
            out = func(*args, **kwargs)

        # Block until the graph materializes and reraise exceptions. This allows
        # `pytest.raises` and `pytest.warns` to work as expected. Note that this would
        # not work on scheduler='distributed', as it would not block.
        arrays, rest = pickle_flatten(out, da.Array)
        arrays = dask.persist(arrays, scheduler="threads")[0]  # type: ignore[attr-defined,no-untyped-call]  # pyright: ignore[reportPrivateImportUsage]
        return pickle_unflatten(arrays, rest)  # pyright: ignore[reportUnknownArgumentType]

    return wrapper


def _require_numpy() -> ModuleType:  # numpydoc ignore=RT01
    """
    Import and return `numpy` if it is available, otherwise raise informative error.
    """
    try:
        import numpy as np
    except ImportError as e:
        msg = (
            "The assertion functions of `xpx.testing` require the numpy module "
            "to be importable in the Python environment."
        )
        raise ImportError(msg) from e

    return np


def _check_ns_shape_dtype(
    actual: Array,
    desired: Array,
    check_dtype: bool,
    check_shape: bool,
    check_scalar: bool,
    xp: ModuleType | None = None,
) -> tuple[Array, Array, ModuleType, ModuleType]:  # numpydoc ignore=RT03
    """
    Assert that namespace, shape and dtype of the two arrays match.

    Parameters
    ----------
    actual : Array
        The array produced by the tested function.
    desired : Array
        The expected array (typically hardcoded).
    check_dtype : bool, default: True
        Whether to check agreement between actual and desired dtypes.
    check_shape : bool, default: True
        Whether to check agreement between actual and desired shapes.
    check_scalar : bool, default: False
        NumPy only: whether to check agreement between actual and desired types -
        0d array vs scalar.
    xp : array_namespace, optional
        A standard-compatible namespace which `actual` and `desired` must match.

    Returns
    -------
    Actual array, desired array, their array namespace, the numpy module.
    """
    np = _require_numpy()

    actual_xp = array_namespace(actual)  # Raises on Python scalars and lists

    if xp is not None:
        _msg = (
            "Namespace of actual array does not match the `xp` argument.\n"
            f"Actual array's namespace: {actual_xp.__name__}\n"
            f"Expected namespace: {xp.__name__}."
        )
        assert actual_xp == xp, _msg
        desired_xp = xp
    else:
        desired_xp = array_namespace(desired)
        _msg = (
            "Namespaces of actual and desired arrays do not match.\n"
            f"Actual: {actual_xp.__name__}\n"
            f"Desired: {desired_xp.__name__}."
        )
        assert actual_xp == desired_xp, _msg

    if is_numpy_namespace(actual_xp) and check_scalar:
        # only NumPy distinguishes between scalars and arrays; we do if check_scalar.
        _msg = (
            "array-ness does not match:\n Actual: "
            f"{type(actual)}\n Desired: {type(desired)}"
        )
        assert np.isscalar(actual) == np.isscalar(desired), _msg

    # Dask uses nan instead of None for unknown shapes
    actual_shape = cast(tuple[float, ...], actual.shape)
    desired_shape = cast(tuple[float, ...], desired.shape)
    assert None not in actual_shape  # Requires explicit support
    assert None not in desired_shape

    if is_dask_namespace(desired_xp):
        if any(math.isnan(i) for i in actual_shape):
            actual.compute_chunk_sizes()  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]
            actual_shape = cast(tuple[float, ...], actual.shape)
        if any(math.isnan(i) for i in desired_shape):
            desired.compute_chunk_sizes()  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]
            desired_shape = cast(tuple[float, ...], desired.shape)

    if check_shape:
        msg = f"shapes do not match: {actual_shape} != {desired_shape}"
        assert actual_shape == desired_shape, msg
    elif desired.ndim > 0:
        # Ignore shape, but check flattened size. This is normally done by
        # np.testing.assert_array_equal etc even when strict=False, but not for
        # non-materializable arrays.
        # This check excludes 0d arrays as they are special-cased in NumPy.
        actual_size = math.prod(actual_shape)
        desired_size = math.prod(desired_shape)
        msg = f"sizes do not match: {actual_size} != {desired_size}"
        assert actual_size == desired_size, msg

    desired = desired_xp.asarray(desired)
    if check_dtype:
        msg = f"dtypes do not match: {actual.dtype} != {desired.dtype}"
        assert actual.dtype == desired.dtype, msg
    desired = desired_xp.broadcast_to(desired, actual_shape)
    return actual, desired, desired_xp, np


def _is_materializable(x: Array) -> bool:  # numpydoc ignore=PR01,RT01
    """
    Return True if you can call `as_numpy_array(x)`; False otherwise.
    """
    # Important: here we assume that we're not tracing -
    # e.g. we're not inside `jax.jit`` nor `cupy.cuda.Stream.begin_capture`.
    return not is_torch_array(x) or x.device.type != "meta"  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]


def _as_numpy_array(  # numpydoc ignore=PR01,RT01
    array: Array, *, xp: ModuleType
) -> np.typing.NDArray[Any]:
    """
    Convert array to NumPy, bypassing GPU-CPU transfer guards and densification guards.
    """
    np = _require_numpy()
    if is_cupy_namespace(xp):
        return xp.asnumpy(array)
    if is_pydata_sparse_namespace(xp):
        return array.todense()  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]

    if is_torch_namespace(xp):
        array = cast(Array, array.resolve_conj())  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]
        array = to_device(array, "cpu")
    if is_array_api_strict_namespace(xp):
        cpu: Device = xp.Device("CPU_DEVICE")
        array = to_device(array, cpu)
    if is_jax_namespace(xp):
        import jax

        # Note: only needed if the transfer guard is enabled
        cpu = cast(Device, jax.devices("cpu")[0])
        array = to_device(array, cpu)

    if hasattr(array, "__dlpack__"):
        try:
            return np.from_dlpack(array)
        except (TypeError, BufferError):
            pass

    return np.asarray(array)


def assert_close(
    actual: Array,
    desired: Array,
    *,
    rtol: float | Array | None = None,
    atol: float | Array = 0,
    equal_nan: bool = True,
    err_msg: str = "",
    verbose: bool = True,
    check_dtype: bool = True,
    check_shape: bool = True,
    check_scalar: bool = False,
    xp: ModuleType | None = None,
) -> None:
    """
    Check that two arrays are close, up to tolerance ``atol + rtol * abs(desired)``.

    This is an interface to :func:`numpy.testing.assert_allclose` which accepts
    any standard-compatible array and performs additional array namespace,
    shape, and dtype checks.

    Parameters
    ----------
    actual : Array
        The array produced by the tested function.
    desired : Array
        The expected array (typically hardcoded).
    rtol : float or Array, optional
        Relative tolerance. Default: dtype-dependent.
    atol : float or Array, optional
        Absolute tolerance. Default: 0.
    equal_nan : bool, default: True
        Whether to consider NaNs in corresponding locations as equal.
    err_msg : str, optional
        Error message to display on failure.
    verbose : bool, default: True
        Whether to include the conflicting arrays in the error message on failure.
    check_dtype : bool, default: True
        Whether to check agreement between actual and desired dtypes.
    check_shape : bool, default: True
        Whether to check agreement between actual and desired shapes.
    check_scalar : bool, default: False
        NumPy only: whether to check agreement between actual and desired types —
        0-D :class:`numpy.ndarray` vs scalar (e.g. :class:`numpy.double`).
    xp : array_namespace, optional
        A standard-compatible namespace which `actual` and `desired` must match.

    Raises
    ------
    AssertionError
        If `actual` and `desired` are not equal up to the defined tolerance.

    ImportError
        If :mod:`numpy` is not importable in the Python environment.

    See Also
    --------
    assert_equal : Similar function for exact equality checks.
    array_api_extra.isclose : Similar function checking closeness, returning a bool.
    numpy.testing.assert_allclose : Similar function for NumPy arrays.

    Notes
    -----
    The default `atol` and `rtol` differ from ``xp.all(xpx.isclose(a, b))``.
    For inexact dtypes, the default `rtol` is
    ``xp.finfo(actual.dtype).eps ** 0.5 * 4``, which for ``float64`` is roughly halfway
    between :math:`\\sqrt{\\epsilon}` and the default for
    :func:`numpy.testing.assert_allclose`, ``1e-7``.
    This gives a more reasonable default for lower precision dtypes,
    for example approximately ``1e-3`` for ``float32``.
    For exact dtypes, the default ``1e-7`` is used.

    Array arguments to `atol` and `rtol` must be valid input to :class:`float`.
    """
    __tracebackhide__ = True
    actual, desired, xp, np = _check_ns_shape_dtype(
        actual, desired, check_dtype, check_shape, check_scalar, xp
    )
    if not _is_materializable(actual):
        return

    if rtol is None:
        if xp.isdtype(actual.dtype, ("real floating", "complex floating")):
            # multiplier of 4 is used as for `np.float64` this puts the default `rtol`
            # roughly half way between sqrt(eps) and the default for
            # `numpy.testing.assert_allclose`, 1e-7
            rtol = xp.finfo(actual.dtype).eps ** 0.5 * 4
        else:
            rtol = 1e-7
    else:
        rtol = float(rtol)

    atol = float(atol)

    actual_np = _as_numpy_array(actual, xp=xp)
    desired_np = _as_numpy_array(desired, xp=xp)
    np.testing.assert_allclose(
        actual_np,
        desired_np,
        rtol=rtol,
        atol=atol,
        equal_nan=equal_nan,
        err_msg=err_msg,
        verbose=verbose,
    )


def assert_equal(
    actual: Array,
    desired: Array,
    *,
    err_msg: str = "",
    verbose: bool = True,
    check_dtype: bool = True,
    check_shape: bool = True,
    check_scalar: bool = False,
    xp: ModuleType | None = None,
) -> None:
    """
    Check that two arrays are equal.

    This is an interface to :func:`numpy.testing.assert_array_equal` which accepts
    any standard-compatible array and performs additional array namespace,
    shape, and dtype checks.

    Parameters
    ----------
    actual : Array
        The array produced by the tested function.
    desired : Array
        The expected array (typically hardcoded).
    err_msg : str, optional
        Error message to display on failure.
    verbose : bool, default: True
        Whether to include the conflicting arrays in the error message on failure.
    check_dtype : bool, default: True
        Whether to check agreement between actual and desired dtypes.
    check_shape : bool, default: True
        Whether to check agreement between actual and desired shapes.
    check_scalar : bool, default: False
        NumPy only: whether to check agreement between actual and desired types —
        0-D :class:`numpy.ndarray` vs scalar (e.g. :class:`numpy.double`).
    xp : array_namespace, optional
        A standard-compatible namespace which `actual` and `desired` must match.

    Raises
    ------
    AssertionError
        If `actual` and `desired` are not equal.

    ImportError
        If :mod:`numpy` is not importable in the Python environment.

    See Also
    --------
    assert_close : Similar function for inexact equality checks.
    numpy.testing.assert_array_equal : Similar function for NumPy arrays.
    """
    __tracebackhide__ = True
    actual, desired, xp, np = _check_ns_shape_dtype(
        actual, desired, check_dtype, check_shape, check_scalar, xp
    )
    if not _is_materializable(actual):
        return
    actual_np = _as_numpy_array(actual, xp=xp)
    desired_np = _as_numpy_array(desired, xp=xp)
    np.testing.assert_array_equal(
        actual_np, desired_np, err_msg=err_msg, verbose=verbose
    )


def assert_less(
    x: Array,
    y: Array,
    *,
    err_msg: str = "",
    verbose: bool = True,
    check_dtype: bool = True,
    check_shape: bool = True,
    check_scalar: bool = False,
    xp: ModuleType | None = None,
) -> None:
    """
    Check that two arrays are ordered by less than.

    This is an interface to :func:`numpy.testing.assert_array_less` which accepts
    any standard-compatible array and performs additional array namespace,
    shape, and dtype checks.

    Parameters
    ----------
    x, y : Array
        Array to compare according to ``x < y`` (elementwise).
    err_msg : str, optional
        Error message to display on failure.
    verbose : bool, default: True
        Whether to include the conflicting arrays in the error message on failure.
    check_dtype : bool, default: True
        Whether to check agreement between the dtypes of `x` and `y`.
    check_shape : bool, default: True
        Whether to check agreement between the shapes of `x` and `y`.
    check_scalar : bool, default: False
        NumPy only: whether to check agreement between actual and desired types —
        0-D :class:`numpy.ndarray` vs scalar (e.g. :class:`numpy.double`).
    xp : array_namespace, optional
        A standard-compatible namespace which `x` and `y` must match.

    Raises
    ------
    AssertionError
        If `x` is not strictly smaller than `y`, elementwise.

    ImportError
        If :mod:`numpy` is not importable in the Python environment.

    See Also
    --------
    assert_close : Similar function for inexact equality checks.
    numpy.testing.assert_array_less : Similar function for NumPy arrays.
    """
    __tracebackhide__ = True
    x, y, xp, np = _check_ns_shape_dtype(
        x, y, check_dtype, check_shape, check_scalar, xp
    )
    if not _is_materializable(x):
        return
    x_np = _as_numpy_array(x, xp=xp)
    y_np = _as_numpy_array(y, xp=xp)
    np.testing.assert_array_less(x_np, y_np, err_msg=err_msg, verbose=verbose)


def assert_close_nulp(
    actual: Array,
    desired: Array,
    *,
    nulp: int = 1,
    check_dtype: bool = True,
    check_shape: bool = True,
    check_scalar: bool = False,
    xp: ModuleType | None = None,
) -> None:
    """
    Compare two arrays relatively to their spacing.

    This is an interface to :func:`numpy.testing.assert_array_almost_equal_nulp`
    which accepts any standard-compatible array and performs
    additional array namespace, shape, and dtype checks.

    Parameters
    ----------
    actual : Array
        The array produced by the tested function.
    desired : Array
        The expected array (typically hardcoded).
    nulp : int, optional
        The maximum number of units in the last place
        for the tolerance check. Default: ``1``.
    check_dtype : bool, default: True
        Whether to check agreement between actual and desired dtypes.
    check_shape : bool, default: True
        Whether to check agreement between actual and desired shapes.
    check_scalar : bool, default: False
        NumPy only: whether to check agreement between actual and desired types —
        0-D :class:`numpy.ndarray` vs scalar (e.g. :class:`numpy.double`).
    xp : array_namespace, optional
        A standard-compatible namespace which `actual` and `desired` must match.

    Raises
    ------
    AssertionError
        If the spacing between `actual` and `desired` for one or more elements is \
        larger than `nulp`.

    ImportError
        If :mod:`numpy` is not importable in the Python environment.

    See Also
    --------
    assert_close : Similar function for inexact equality checks.
    numpy.spacing : Spacing calculation for NumPy arrays.
    numpy.testing.assert_array_almost_equal_nulp : Similar function for NumPy arrays.

    Notes
    -----
    This is a relatively robust method to compare two arrays whose amplitude is
    variable.

    An assertion is raised if the following condition is not met::

        abs(actual - desired) <= nulp * spacing(maximum(abs(actual), abs(desired)))

    where ``spacing(x)`` is the distance between ``x`` and the nearest adjacent number
    representable by in the data type of ``x``.
    """
    actual, desired, xp, np = _check_ns_shape_dtype(
        actual, desired, check_dtype, check_shape, check_scalar, xp
    )
    if not _is_materializable(actual):
        return
    actual_np = _as_numpy_array(actual, xp=xp)
    desired_np = _as_numpy_array(desired, xp=xp)
    np.testing.assert_array_almost_equal_nulp(actual_np, desired_np, nulp=nulp)
