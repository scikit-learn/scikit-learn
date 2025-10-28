from __future__ import annotations

from warnings import warn


def find_stacklevel() -> int:
    """Find the first place in the stack that is not inside narwhals.

    Taken from:
    https://github.com/pandas-dev/pandas/blob/ab89c53f48df67709a533b6a95ce3d911871a0a8/pandas/util/_exceptions.py#L30-L51
    """
    import inspect
    from pathlib import Path

    import narwhals as nw

    pkg_dir = str(Path(nw.__file__).parent)

    # https://stackoverflow.com/questions/17407119/python-inspect-stack-is-slow
    frame = inspect.currentframe()
    n = 0
    try:
        while frame:
            fname = inspect.getfile(frame)
            if fname.startswith(pkg_dir) or (
                (qualname := getattr(frame.f_code, "co_qualname", None))
                # ignore @singledispatch wrappers
                and qualname.startswith("singledispatch.")
            ):
                frame = frame.f_back
                n += 1
            else:  # pragma: no cover
                break
        else:  # pragma: no cover
            pass
    finally:
        # https://docs.python.org/3/library/inspect.html
        # > Though the cycle detector will catch these, destruction of the frames
        # > (and local variables) can be made deterministic by removing the cycle
        # > in a finally clause.
        del frame
    return n


def issue_deprecation_warning(message: str, _version: str) -> None:  # pragma: no cover
    """Issue a deprecation warning.

    Arguments:
        message: The message associated with the warning.
        _version: Narwhals version when the warning was introduced. Just used for internal
            bookkeeping.
    """
    warn(message=message, category=DeprecationWarning, stacklevel=find_stacklevel())


def issue_warning(message: str, category: type[Warning]) -> None:
    warn(message=message, category=category, stacklevel=find_stacklevel())
