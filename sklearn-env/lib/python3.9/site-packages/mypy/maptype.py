from typing import Dict, List

from mypy.expandtype import expand_type
from mypy.nodes import TypeInfo
from mypy.types import Type, TypeVarId, Instance, AnyType, TypeOfAny, ProperType


def map_instance_to_supertype(instance: Instance,
                              superclass: TypeInfo) -> Instance:
    """Produce a supertype of `instance` that is an Instance
    of `superclass`, mapping type arguments up the chain of bases.

    If `superclass` is not a nominal superclass of `instance.type`,
    then all type arguments are mapped to 'Any'.
    """
    if instance.type == superclass:
        # Fast path: `instance` already belongs to `superclass`.
        return instance

    if not superclass.type_vars:
        # Fast path: `superclass` has no type variables to map to.
        return Instance(superclass, [])

    return map_instance_to_supertypes(instance, superclass)[0]


def map_instance_to_supertypes(instance: Instance,
                               supertype: TypeInfo) -> List[Instance]:
    # FIX: Currently we should only have one supertype per interface, so no
    #      need to return an array
    result: List[Instance] = []
    for path in class_derivation_paths(instance.type, supertype):
        types = [instance]
        for sup in path:
            a: List[Instance] = []
            for t in types:
                a.extend(map_instance_to_direct_supertypes(t, sup))
            types = a
        result.extend(types)
    if result:
        return result
    else:
        # Nothing. Presumably due to an error. Construct a dummy using Any.
        any_type = AnyType(TypeOfAny.from_error)
        return [Instance(supertype, [any_type] * len(supertype.type_vars))]


def class_derivation_paths(typ: TypeInfo,
                           supertype: TypeInfo) -> List[List[TypeInfo]]:
    """Return an array of non-empty paths of direct base classes from
    type to supertype.  Return [] if no such path could be found.

      InterfaceImplementationPaths(A, B) == [[B]] if A inherits B
      InterfaceImplementationPaths(A, C) == [[B, C]] if A inherits B and
                                                        B inherits C
    """
    # FIX: Currently we might only ever have a single path, so this could be
    #      simplified
    result: List[List[TypeInfo]] = []

    for base in typ.bases:
        btype = base.type
        if btype == supertype:
            result.append([btype])
        else:
            # Try constructing a longer path via the base class.
            for path in class_derivation_paths(btype, supertype):
                result.append([btype] + path)

    return result


def map_instance_to_direct_supertypes(instance: Instance,
                                      supertype: TypeInfo) -> List[Instance]:
    # FIX: There should only be one supertypes, always.
    typ = instance.type
    result: List[Instance] = []

    for b in typ.bases:
        if b.type == supertype:
            env = instance_to_type_environment(instance)
            t = expand_type(b, env)
            assert isinstance(t, ProperType)
            assert isinstance(t, Instance)
            result.append(t)

    if result:
        return result
    else:
        # Relationship with the supertype not specified explicitly. Use dynamic
        # type arguments implicitly.
        any_type = AnyType(TypeOfAny.unannotated)
        return [Instance(supertype, [any_type] * len(supertype.type_vars))]


def instance_to_type_environment(instance: Instance) -> Dict[TypeVarId, Type]:
    """Given an Instance, produce the resulting type environment for type
    variables bound by the Instance's class definition.

    An Instance is a type application of a class (a TypeInfo) to its
    required number of type arguments.  So this environment consists
    of the class's type variables mapped to the Instance's actual
    arguments.  The type variables are mapped by their `id`.

    """
    return {binder.id: arg for binder, arg in zip(instance.type.defn.type_vars, instance.args)}
