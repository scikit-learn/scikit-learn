"""Maintain a mapping from mypy concepts to IR/compiled concepts."""

from typing import Dict, Optional

from mypy.nodes import FuncDef, TypeInfo, SymbolNode, ArgKind, ARG_STAR, ARG_STAR2
from mypy.types import (
    Instance, Type, CallableType, LiteralType, TypedDictType, UnboundType, PartialType,
    UninhabitedType, Overloaded, UnionType, TypeType, AnyType, NoneTyp, TupleType, TypeVarType,
    get_proper_type
)

from mypyc.ir.rtypes import (
    RType, RUnion, RTuple, RInstance, object_rprimitive, dict_rprimitive, tuple_rprimitive,
    none_rprimitive, int_rprimitive, float_rprimitive, str_rprimitive, bool_rprimitive,
    list_rprimitive, set_rprimitive, range_rprimitive, bytes_rprimitive
)
from mypyc.ir.func_ir import FuncSignature, FuncDecl, RuntimeArg
from mypyc.ir.class_ir import ClassIR


class Mapper:
    """Keep track of mappings from mypy concepts to IR concepts.

    For example, we keep track of how the mypy TypeInfos of compiled
    classes map to class IR objects.

    This state is shared across all modules being compiled in all
    compilation groups.
    """

    def __init__(self, group_map: Dict[str, Optional[str]]) -> None:
        self.group_map = group_map
        self.type_to_ir: Dict[TypeInfo, ClassIR] = {}
        self.func_to_decl: Dict[SymbolNode, FuncDecl] = {}

    def type_to_rtype(self, typ: Optional[Type]) -> RType:
        if typ is None:
            return object_rprimitive

        typ = get_proper_type(typ)
        if isinstance(typ, Instance):
            if typ.type.fullname == 'builtins.int':
                return int_rprimitive
            elif typ.type.fullname == 'builtins.float':
                return float_rprimitive
            elif typ.type.fullname == 'builtins.bool':
                return bool_rprimitive
            elif typ.type.fullname == 'builtins.str':
                return str_rprimitive
            elif typ.type.fullname == 'builtins.bytes':
                return bytes_rprimitive
            elif typ.type.fullname == 'builtins.list':
                return list_rprimitive
            # Dict subclasses are at least somewhat common and we
            # specifically support them, so make sure that dict operations
            # get optimized on them.
            elif any(cls.fullname == 'builtins.dict' for cls in typ.type.mro):
                return dict_rprimitive
            elif typ.type.fullname == 'builtins.set':
                return set_rprimitive
            elif typ.type.fullname == 'builtins.tuple':
                return tuple_rprimitive  # Varying-length tuple
            elif typ.type.fullname == 'builtins.range':
                return range_rprimitive
            elif typ.type in self.type_to_ir:
                inst = RInstance(self.type_to_ir[typ.type])
                # Treat protocols as Union[protocol, object], so that we can do fast
                # method calls in the cases where the protocol is explicitly inherited from
                # and fall back to generic operations when it isn't.
                if typ.type.is_protocol:
                    return RUnion([inst, object_rprimitive])
                else:
                    return inst
            else:
                return object_rprimitive
        elif isinstance(typ, TupleType):
            # Use our unboxed tuples for raw tuples but fall back to
            # being boxed for NamedTuple.
            if typ.partial_fallback.type.fullname == 'builtins.tuple':
                return RTuple([self.type_to_rtype(t) for t in typ.items])
            else:
                return tuple_rprimitive
        elif isinstance(typ, CallableType):
            return object_rprimitive
        elif isinstance(typ, NoneTyp):
            return none_rprimitive
        elif isinstance(typ, UnionType):
            return RUnion([self.type_to_rtype(item)
                           for item in typ.items])
        elif isinstance(typ, AnyType):
            return object_rprimitive
        elif isinstance(typ, TypeType):
            return object_rprimitive
        elif isinstance(typ, TypeVarType):
            # Erase type variable to upper bound.
            # TODO: Erase to union if object has value restriction?
            return self.type_to_rtype(typ.upper_bound)
        elif isinstance(typ, PartialType):
            assert typ.var.type is not None
            return self.type_to_rtype(typ.var.type)
        elif isinstance(typ, Overloaded):
            return object_rprimitive
        elif isinstance(typ, TypedDictType):
            return dict_rprimitive
        elif isinstance(typ, LiteralType):
            return self.type_to_rtype(typ.fallback)
        elif isinstance(typ, (UninhabitedType, UnboundType)):
            # Sure, whatever!
            return object_rprimitive

        # I think we've covered everything that is supposed to
        # actually show up, so anything else is a bug somewhere.
        assert False, 'unexpected type %s' % type(typ)

    def get_arg_rtype(self, typ: Type, kind: ArgKind) -> RType:
        if kind == ARG_STAR:
            return tuple_rprimitive
        elif kind == ARG_STAR2:
            return dict_rprimitive
        else:
            return self.type_to_rtype(typ)

    def fdef_to_sig(self, fdef: FuncDef) -> FuncSignature:
        if isinstance(fdef.type, CallableType):
            arg_types = [self.get_arg_rtype(typ, kind)
                         for typ, kind in zip(fdef.type.arg_types, fdef.type.arg_kinds)]
            arg_pos_onlys = [name is None for name in fdef.type.arg_names]
            ret = self.type_to_rtype(fdef.type.ret_type)
        else:
            # Handle unannotated functions
            arg_types = [object_rprimitive for arg in fdef.arguments]
            arg_pos_onlys = [arg.pos_only for arg in fdef.arguments]
            # We at least know the return type for __init__ methods will be None.
            is_init_method = fdef.name == '__init__' and bool(fdef.info)
            if is_init_method:
                ret = none_rprimitive
            else:
                ret = object_rprimitive

        # mypyc FuncSignatures (unlike mypy types) want to have a name
        # present even when the argument is position only, since it is
        # the sole way that FuncDecl arguments are tracked. This is
        # generally fine except in some cases (like for computing
        # init_sig) we need to produce FuncSignatures from a
        # deserialized FuncDef that lacks arguments. We won't ever
        # need to use those inside of a FuncIR, so we just make up
        # some crap.
        if hasattr(fdef, 'arguments'):
            arg_names = [arg.variable.name for arg in fdef.arguments]
        else:
            arg_names = [name or '' for name in fdef.arg_names]

        args = [RuntimeArg(arg_name, arg_type, arg_kind, arg_pos_only)
                for arg_name, arg_kind, arg_type, arg_pos_only
                in zip(arg_names, fdef.arg_kinds, arg_types, arg_pos_onlys)]

        # We force certain dunder methods to return objects to support letting them
        # return NotImplemented. It also avoids some pointless boxing and unboxing,
        # since tp_richcompare needs an object anyways.
        if fdef.name in ('__eq__', '__ne__', '__lt__', '__gt__', '__le__', '__ge__'):
            ret = object_rprimitive
        return FuncSignature(args, ret)
