from __future__ import annotations

from mypyc.ir.ops import GetElementPtr, LoadMem, Value
from mypyc.ir.rtypes import PyVarObject, c_pyssize_t_rprimitive
from mypyc.irbuild.ll_builder import LowLevelIRBuilder
from mypyc.lower.registry import lower_primitive_op


@lower_primitive_op("var_object_size")
def var_object_size(builder: LowLevelIRBuilder, args: list[Value], line: int) -> Value:
    elem_address = builder.add(GetElementPtr(args[0], PyVarObject, "ob_size"))
    return builder.add(LoadMem(c_pyssize_t_rprimitive, elem_address))
