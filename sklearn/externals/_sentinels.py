# type: ignore
"""
Copied from https://github.com/taleinat/python-stdlib-sentinels
PEP-0661: Status: Draft
"""
import sys as _sys
from typing import Optional


__all__ = ["sentinel"]


def sentinel(
    name: str,
    repr: Optional[str] = None,
    module: Optional[str] = None,
):
    """Create a unique sentinel object.

    *name* should be the fully-qualified name of the variable to which the
    return value shall be assigned.

    *repr*, if supplied, will be used for the repr of the sentinel object.
    If not provided, "<name>" will be used (with any leading class names
    removed).

    *module*, if supplied, will be used as the module name for the purpose
    of setting a unique name for the sentinels unique class.  The class is
    set as an attribute of this name on the "sentinels" module, so that it
    may be found by the pickling mechanism.  In most cases, the module name
    does not need to be provided, and it will be found by inspecting the
    stack frame.
    """
    name = _sys.intern(str(name))
    repr = repr or f'<{name.rsplit(".", 1)[-1]}>'

    if module is None:
        try:
            module = _get_parent_frame().f_globals.get("__name__", "__main__")
        except (AttributeError, ValueError):
            pass
    class_name = _sys.intern(_get_class_name(name, module))

    class_namespace = {
        "__repr__": lambda self: repr,
    }
    cls = type(class_name, (), class_namespace)
    cls.__module__ = __name__
    globals()[class_name] = cls

    sentinel = cls()

    def __new__(cls):
        return sentinel

    __new__.__qualname__ = f"{class_name}.__new__"
    cls.__new__ = __new__

    return sentinel


if hasattr(_sys, "_getframe"):
    _get_parent_frame = lambda: _sys._getframe(2)
else:  # pragma: no cover

    def _get_parent_frame():
        """Return the frame object for the caller's parent stack frame."""
        try:
            raise Exception
        except Exception:
            return _sys.exc_info()[2].tb_frame.f_back.f_back


def _get_class_name(
    sentinel_qualname: str,
    module_name: Optional[str] = None,
) -> str:
    return (
        "_sentinel_type__"
        f'{module_name.replace(".", "_") + "__" if module_name else ""}'
        f'{sentinel_qualname.replace(".", "_")}'
    )
