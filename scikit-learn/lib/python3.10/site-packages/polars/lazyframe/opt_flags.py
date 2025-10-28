from __future__ import annotations

import contextlib

from polars._utils.deprecation import issue_deprecation_warning

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PyOptFlags

import inspect
from functools import wraps
from typing import TYPE_CHECKING, Callable, TypeVar

if TYPE_CHECKING:
    import sys

    if sys.version_info >= (3, 10):
        from typing import ParamSpec
    else:
        from typing_extensions import ParamSpec

    P = ParamSpec("P")
    T = TypeVar("T")


class QueryOptFlags:
    """
    The set of the optimizations considered during query optimization.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.
    """

    def __init__(
        self,
        *,
        predicate_pushdown: None | bool = None,
        projection_pushdown: None | bool = None,
        simplify_expression: None | bool = None,
        slice_pushdown: None | bool = None,
        comm_subplan_elim: None | bool = None,
        comm_subexpr_elim: None | bool = None,
        cluster_with_columns: None | bool = None,
        collapse_joins: None | bool = None,
        check_order_observe: None | bool = None,
        fast_projection: None | bool = None,
    ) -> None:
        self._pyoptflags = PyOptFlags.default()
        self.update(
            predicate_pushdown=predicate_pushdown,
            projection_pushdown=projection_pushdown,
            simplify_expression=simplify_expression,
            slice_pushdown=slice_pushdown,
            comm_subplan_elim=comm_subplan_elim,
            comm_subexpr_elim=comm_subexpr_elim,
            cluster_with_columns=cluster_with_columns,
            collapse_joins=collapse_joins,
            check_order_observe=check_order_observe,
            fast_projection=fast_projection,
        )

    @classmethod
    def _from_pyoptflags(self, pyoptflags: PyOptFlags) -> QueryOptFlags:
        optflags = self.__new__(self)
        optflags._pyoptflags = pyoptflags
        return optflags

    @staticmethod
    def none(
        *,
        predicate_pushdown: None | bool = None,
        projection_pushdown: None | bool = None,
        simplify_expression: None | bool = None,
        slice_pushdown: None | bool = None,
        comm_subplan_elim: None | bool = None,
        comm_subexpr_elim: None | bool = None,
        cluster_with_columns: None | bool = None,
        collapse_joins: None | bool = None,
        check_order_observe: None | bool = None,
        fast_projection: None | bool = None,
    ) -> QueryOptFlags:
        """Create new empty set off optimizations."""
        optflags = QueryOptFlags()
        optflags.no_optimizations()
        return optflags.update(
            predicate_pushdown=predicate_pushdown,
            projection_pushdown=projection_pushdown,
            simplify_expression=simplify_expression,
            slice_pushdown=slice_pushdown,
            comm_subplan_elim=comm_subplan_elim,
            comm_subexpr_elim=comm_subexpr_elim,
            cluster_with_columns=cluster_with_columns,
            collapse_joins=collapse_joins,
            check_order_observe=check_order_observe,
            fast_projection=fast_projection,
        )

    def update(
        self,
        *,
        predicate_pushdown: None | bool = None,
        projection_pushdown: None | bool = None,
        simplify_expression: None | bool = None,
        slice_pushdown: None | bool = None,
        comm_subplan_elim: None | bool = None,
        comm_subexpr_elim: None | bool = None,
        cluster_with_columns: None | bool = None,
        collapse_joins: None | bool = None,
        check_order_observe: None | bool = None,
        fast_projection: None | bool = None,
    ) -> QueryOptFlags:
        """Update the current optimization flags."""
        if predicate_pushdown is not None:
            self.predicate_pushdown = predicate_pushdown
        if projection_pushdown is not None:
            self.projection_pushdown = projection_pushdown
        if simplify_expression is not None:
            self.simplify_expression = simplify_expression
        if slice_pushdown is not None:
            self.slice_pushdown = slice_pushdown
        if comm_subplan_elim is not None:
            self.comm_subplan_elim = comm_subplan_elim
        if comm_subexpr_elim is not None:
            self.comm_subexpr_elim = comm_subexpr_elim
        if cluster_with_columns is not None:
            self.cluster_with_columns = cluster_with_columns
        if collapse_joins is not None:
            issue_deprecation_warning(
                "the `collapse_joins` parameter for `QueryOptFlags` is deprecated. "
                "Use `predicate_pushdown` instead.",
                version="1.33.1",
            )
            if not collapse_joins:
                self.predicate_pushdown = False
        if check_order_observe is not None:
            self.check_order_observe = check_order_observe
        if fast_projection is not None:
            self.fast_projection = fast_projection

        return self

    @staticmethod
    def _eager() -> QueryOptFlags:
        """Create new empty set off optimizations."""
        optflags = QueryOptFlags()
        optflags.no_optimizations()
        optflags._pyoptflags.eager = True
        optflags.simplify_expression = True
        return optflags

    def __copy__(self) -> QueryOptFlags:
        return QueryOptFlags._from_pyoptflags(self._pyoptflags.copy())

    def __deepcopy__(self) -> QueryOptFlags:
        return QueryOptFlags._from_pyoptflags(self._pyoptflags.copy())

    def no_optimizations(self) -> None:
        """Remove selected optimizations."""
        self._pyoptflags.no_optimizations()

    @property
    def projection_pushdown(self) -> bool:
        """Only read columns that are used later in the query."""
        return self._pyoptflags.projection_pushdown

    @projection_pushdown.setter
    def projection_pushdown(self, value: bool) -> None:
        self._pyoptflags.projection_pushdown = value

    @property
    def predicate_pushdown(self) -> bool:
        """Apply predicates/filters as early as possible."""
        return self._pyoptflags.predicate_pushdown

    @predicate_pushdown.setter
    def predicate_pushdown(self, value: bool) -> None:
        self._pyoptflags.predicate_pushdown = value

    @property
    def cluster_with_columns(self) -> bool:
        """Cluster sequential `with_columns` calls to independent calls."""
        return self._pyoptflags.cluster_with_columns

    @cluster_with_columns.setter
    def cluster_with_columns(self, value: bool) -> None:
        self._pyoptflags.cluster_with_columns = value

    @property
    def simplify_expression(self) -> bool:
        """Run many expression optimization rules until fixed point."""
        return self._pyoptflags.simplify_expression

    @simplify_expression.setter
    def simplify_expression(self, value: bool) -> None:
        self._pyoptflags.simplify_expression = value

    @property
    def slice_pushdown(self) -> bool:
        """Pushdown slices/limits."""
        return self._pyoptflags.slice_pushdown

    @slice_pushdown.setter
    def slice_pushdown(self, value: bool) -> None:
        self._pyoptflags.slice_pushdown = value

    @property
    def comm_subplan_elim(self) -> bool:
        """Elide duplicate plans and caches their outputs."""
        return self._pyoptflags.comm_subplan_elim

    @comm_subplan_elim.setter
    def comm_subplan_elim(self, value: bool) -> None:
        self._pyoptflags.comm_subplan_elim = value

    @property
    def comm_subexpr_elim(self) -> bool:
        """Elide duplicate expressions and caches their outputs."""
        return self._pyoptflags.comm_subexpr_elim

    @comm_subexpr_elim.setter
    def comm_subexpr_elim(self, value: bool) -> None:
        self._pyoptflags.comm_subexpr_elim = value

    @property
    def check_order_observe(self) -> bool:
        """Do not maintain order if the order would not be observed."""
        return self._pyoptflags.check_order_observe

    @check_order_observe.setter
    def check_order_observe(self, value: bool) -> None:
        self._pyoptflags.check_order_observe = value

    @property
    def fast_projection(self) -> bool:
        """Replace simple projections with a faster inlined projection that skips the expression engine."""  # noqa: W505
        return self._pyoptflags.fast_projection

    @fast_projection.setter
    def fast_projection(self, value: bool) -> None:
        self._pyoptflags.fast_projection = value

    def __str__(self) -> str:
        return f"""
QueryOptFlags {{
    type_coercion: {self._pyoptflags.type_coercion}
    type_check: {self._pyoptflags.type_check}

    predicate_pushdown: {self.predicate_pushdown}
    projection_pushdown: {self.projection_pushdown}
    simplify_expression: {self.simplify_expression}
    slice_pushdown: {self.slice_pushdown}
    comm_subplan_elim: {self.comm_subplan_elim}
    comm_subexpr_elim: {self.comm_subexpr_elim}
    cluster_with_columns: {self.cluster_with_columns}
    check_order_observe: {self.check_order_observe}
    fast_projection: {self.fast_projection}

    eager: {self._pyoptflags.eager}
    streaming: {self._pyoptflags.streaming}
}}
        """.strip()


DEFAULT_QUERY_OPT_FLAGS: QueryOptFlags
try:  # Module not available when building docs
    DEFAULT_QUERY_OPT_FLAGS = QueryOptFlags()
except (ImportError, NameError) as _:
    DEFAULT_QUERY_OPT_FLAGS = ()  # type: ignore[assignment]


def forward_old_opt_flags() -> Callable[[Callable[P, T]], Callable[P, T]]:
    """Decorator to mark to forward the old optimization flags."""

    def helper(f: QueryOptFlags, field_name: str, value: bool) -> QueryOptFlags:  # noqa: FBT001
        setattr(f, field_name, value)
        return f

    def helper_hidden(f: QueryOptFlags, field_name: str, value: bool) -> QueryOptFlags:  # noqa: FBT001
        setattr(f._pyoptflags, field_name, value)
        return f

    def clear_optimizations(f: QueryOptFlags, value: bool) -> QueryOptFlags:  # noqa: FBT001
        if value:
            return QueryOptFlags.none()
        else:
            return f

    def eager(f: QueryOptFlags, value: bool) -> QueryOptFlags:  # noqa: FBT001
        if value:
            return QueryOptFlags._eager()
        else:
            return f

    OLD_OPT_PARAMETERS_MAPPING = {
        "no_optimization": lambda f, v: clear_optimizations(f, v),
        "_eager": lambda f, v: eager(f, v),
        "type_coercion": lambda f, v: helper_hidden(f, "type_coercion", v),
        "_type_check": lambda f, v: helper_hidden(f, "type_check", v),
        "predicate_pushdown": lambda f, v: helper(f, "predicate_pushdown", v),
        "projection_pushdown": lambda f, v: helper(f, "projection_pushdown", v),
        "simplify_expression": lambda f, v: helper(f, "simplify_expression", v),
        "slice_pushdown": lambda f, v: helper(f, "slice_pushdown", v),
        "comm_subplan_elim": lambda f, v: helper(f, "comm_subplan_elim", v),
        "comm_subexpr_elim": lambda f, v: helper(f, "comm_subexpr_elim", v),
        "cluster_with_columns": lambda f, v: helper(f, "cluster_with_columns", v),
        "collapse_joins": lambda f, v: helper(f, "collapse_joins", v),
        "_check_order": lambda f, v: helper(f, "check_order_observe", v),
    }

    def decorate(function: Callable[P, T]) -> Callable[P, T]:
        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            optflags: QueryOptFlags = kwargs.get(
                "optimizations", DEFAULT_QUERY_OPT_FLAGS
            )  # type: ignore[assignment]
            optflags = optflags.__copy__()
            for key in list(kwargs.keys()):
                cb = OLD_OPT_PARAMETERS_MAPPING.get(key)
                if cb is not None:
                    from polars._utils.various import issue_warning

                    message = f"optimization flag `{key}` is deprecated. Please use `optimizations` parameter\n(Deprecated in version 1.30.0)"
                    issue_warning(message, DeprecationWarning)
                    optflags = cb(optflags, kwargs.pop(key))  # type: ignore[no-untyped-call,unused-ignore]

            kwargs["optimizations"] = optflags
            return function(*args, **kwargs)

        wrapper.__signature__ = inspect.signature(function)  # type: ignore[attr-defined]
        return wrapper

    return decorate
