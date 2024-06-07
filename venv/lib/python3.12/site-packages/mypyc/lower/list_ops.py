from __future__ import annotations

from mypyc.common import PLATFORM_SIZE
from mypyc.ir.ops import GetElementPtr, Integer, IntOp, LoadMem, SetMem, Value
from mypyc.ir.rtypes import (
    PyListObject,
    c_pyssize_t_rprimitive,
    object_rprimitive,
    pointer_rprimitive,
)
from mypyc.irbuild.ll_builder import LowLevelIRBuilder
from mypyc.lower.registry import lower_primitive_op


@lower_primitive_op("buf_init_item")
def buf_init_item(builder: LowLevelIRBuilder, args: list[Value], line: int) -> Value:
    """Initialize an item in a buffer of "PyObject *" values at given index.

    This can be used to initialize the data buffer of a freshly allocated list
    object.
    """
    base = args[0]
    index_value = args[1]
    value = args[2]
    assert isinstance(index_value, Integer)
    index = index_value.numeric_value()
    if index == 0:
        ptr = base
    else:
        ptr = builder.add(
            IntOp(
                pointer_rprimitive,
                base,
                Integer(index * PLATFORM_SIZE, c_pyssize_t_rprimitive),
                IntOp.ADD,
                line,
            )
        )
    return builder.add(SetMem(object_rprimitive, ptr, value, line))


@lower_primitive_op("list_items")
def list_items(builder: LowLevelIRBuilder, args: list[Value], line: int) -> Value:
    ob_item_ptr = builder.add(GetElementPtr(args[0], PyListObject, "ob_item", line))
    return builder.add(LoadMem(pointer_rprimitive, ob_item_ptr, line))
