"""Functions to check that serialization round-tripped properly."""

# This file is named test_serialization.py even though it doesn't
# contain its own tests so that pytest will rewrite the asserts...

from typing import Any, Dict, Tuple
from mypy.backports import OrderedDict
from collections.abc import Iterable

from mypyc.ir.ops import DeserMaps
from mypyc.ir.rtypes import RType
from mypyc.ir.func_ir import FuncDecl, FuncIR, FuncSignature
from mypyc.ir.class_ir import ClassIR
from mypyc.ir.module_ir import ModuleIR, deserialize_modules
from mypyc.sametype import is_same_type, is_same_signature


def get_dict(x: Any) -> Dict[str, Any]:
    if hasattr(x, '__mypyc_attrs__'):
        return {k: getattr(x, k) for k in x.__mypyc_attrs__ if hasattr(x, k)}
    else:
        return dict(x.__dict__)


def get_function_dict(x: FuncIR) -> Dict[str, Any]:
    """Get a dict of function attributes safe to compare across serialization"""
    d = get_dict(x)
    d.pop('blocks', None)
    d.pop('env', None)
    return d


def assert_blobs_same(x: Any, y: Any, trail: Tuple[Any, ...]) -> None:
    """Compare two blobs of IR as best we can.

    FuncDecls, FuncIRs, and ClassIRs are compared by fullname to avoid
    infinite recursion.
    (More detailed comparisons should be done manually.)

    Types and signatures are compared using mypyc.sametype.

    Containers are compared recursively.

    Anything else is compared with ==.

    The `trail` argument is used in error messages.
    """

    assert type(x) is type(y), ("Type mismatch at {}".format(trail), type(x), type(y))
    if isinstance(x, (FuncDecl, FuncIR, ClassIR)):
        assert x.fullname == y.fullname, "Name mismatch at {}".format(trail)
    elif isinstance(x, OrderedDict):
        assert len(x.keys()) == len(y.keys()), "Keys mismatch at {}".format(trail)
        for (xk, xv), (yk, yv) in zip(x.items(), y.items()):
            assert_blobs_same(xk, yk, trail + ("keys",))
            assert_blobs_same(xv, yv, trail + (xk,))
    elif isinstance(x, dict):
        assert x.keys() == y.keys(), "Keys mismatch at {}".format(trail)
        for k in x.keys():
            assert_blobs_same(x[k], y[k], trail + (k,))
    elif isinstance(x, Iterable) and not isinstance(x, str):
        for i, (xv, yv) in enumerate(zip(x, y)):
            assert_blobs_same(xv, yv, trail + (i,))
    elif isinstance(x, RType):
        assert is_same_type(x, y), "RType mismatch at {}".format(trail)
    elif isinstance(x, FuncSignature):
        assert is_same_signature(x, y), "Signature mismatch at {}".format(trail)
    else:
        assert x == y, "Value mismatch at {}".format(trail)


def assert_modules_same(ir1: ModuleIR, ir2: ModuleIR) -> None:
    """Assert that two module IRs are the same (*).

    * Or rather, as much as we care about preserving across
    serialization.  We drop the actual IR bodies of functions but try
    to preserve everything else.
    """
    assert ir1.fullname == ir2.fullname

    assert ir1.imports == ir2.imports

    for cls1, cls2 in zip(ir1.classes, ir2.classes):
        assert_blobs_same(get_dict(cls1), get_dict(cls2), (ir1.fullname, cls1.fullname))

    for fn1, fn2 in zip(ir1.functions, ir2.functions):
        assert_blobs_same(get_function_dict(fn1), get_function_dict(fn2),
                          (ir1.fullname, fn1.fullname))
        assert_blobs_same(get_dict(fn1.decl), get_dict(fn2.decl),
                          (ir1.fullname, fn1.fullname))

    assert_blobs_same(ir1.final_names, ir2.final_names, (ir1.fullname, 'final_names'))


def check_serialization_roundtrip(irs: Dict[str, ModuleIR]) -> None:
    """Check that we can serialize modules out and deserialize them to the same thing."""
    serialized = {k: ir.serialize() for k, ir in irs.items()}

    ctx = DeserMaps({}, {})
    irs2 = deserialize_modules(serialized, ctx)
    assert irs.keys() == irs2.keys()

    for k in irs:
        assert_modules_same(irs[k], irs2[k])
