"""Type inference constraints."""

from typing import TYPE_CHECKING, Iterable, List, Optional, Sequence
from typing_extensions import Final

from mypy.types import (
    CallableType, Type, TypeVisitor, UnboundType, AnyType, NoneType, TypeVarType, Instance,
    TupleType, TypedDictType, UnionType, Overloaded, ErasedType, PartialType, DeletedType,
    UninhabitedType, TypeType, TypeVarId, TypeQuery, is_named_instance, TypeOfAny, LiteralType,
    ProperType, ParamSpecType, get_proper_type, TypeAliasType, is_union_with_any,
    callable_with_ellipsis
)
from mypy.maptype import map_instance_to_supertype
import mypy.subtypes
import mypy.sametypes
import mypy.typeops
from mypy.erasetype import erase_typevars
from mypy.nodes import COVARIANT, CONTRAVARIANT, ArgKind
from mypy.argmap import ArgTypeExpander
from mypy.typestate import TypeState

if TYPE_CHECKING:
    from mypy.infer import ArgumentInferContext

SUBTYPE_OF: Final = 0
SUPERTYPE_OF: Final = 1


class Constraint:
    """A representation of a type constraint.

    It can be either T <: type or T :> type (T is a type variable).
    """

    type_var: TypeVarId
    op = 0           # SUBTYPE_OF or SUPERTYPE_OF
    target: Type

    def __init__(self, type_var: TypeVarId, op: int, target: Type) -> None:
        self.type_var = type_var
        self.op = op
        self.target = target

    def __repr__(self) -> str:
        op_str = '<:'
        if self.op == SUPERTYPE_OF:
            op_str = ':>'
        return '{} {} {}'.format(self.type_var, op_str, self.target)


def infer_constraints_for_callable(
        callee: CallableType,
        arg_types: Sequence[Optional[Type]],
        arg_kinds: List[ArgKind],
        formal_to_actual: List[List[int]],
        context: 'ArgumentInferContext') -> List[Constraint]:
    """Infer type variable constraints for a callable and actual arguments.

    Return a list of constraints.
    """
    constraints: List[Constraint] = []
    mapper = ArgTypeExpander(context)

    for i, actuals in enumerate(formal_to_actual):
        for actual in actuals:
            actual_arg_type = arg_types[actual]
            if actual_arg_type is None:
                continue

            actual_type = mapper.expand_actual_type(actual_arg_type, arg_kinds[actual],
                                                    callee.arg_names[i], callee.arg_kinds[i])
            c = infer_constraints(callee.arg_types[i], actual_type, SUPERTYPE_OF)
            constraints.extend(c)

    return constraints


def infer_constraints(template: Type, actual: Type,
                      direction: int) -> List[Constraint]:
    """Infer type constraints.

    Match a template type, which may contain type variable references,
    recursively against a type which does not contain (the same) type
    variable references. The result is a list of type constrains of
    form 'T is a supertype/subtype of x', where T is a type variable
    present in the template and x is a type without reference to type
    variables present in the template.

    Assume T and S are type variables. Now the following results can be
    calculated (read as '(template, actual) --> result'):

      (T, X)            -->  T :> X
      (X[T], X[Y])      -->  T <: Y and T :> Y
      ((T, T), (X, Y))  -->  T :> X and T :> Y
      ((T, S), (X, Y))  -->  T :> X and S :> Y
      (X[T], Any)       -->  T <: Any and T :> Any

    The constraints are represented as Constraint objects.
    """
    if any(get_proper_type(template) == get_proper_type(t) for t in TypeState._inferring):
        return []
    if isinstance(template, TypeAliasType) and template.is_recursive:
        # This case requires special care because it may cause infinite recursion.
        TypeState._inferring.append(template)
        res = _infer_constraints(template, actual, direction)
        TypeState._inferring.pop()
        return res
    return _infer_constraints(template, actual, direction)


def _infer_constraints(template: Type, actual: Type,
                       direction: int) -> List[Constraint]:

    orig_template = template
    template = get_proper_type(template)
    actual = get_proper_type(actual)

    # Type inference shouldn't be affected by whether union types have been simplified.
    # We however keep any ErasedType items, so that the caller will see it when using
    # checkexpr.has_erased_component().
    if isinstance(template, UnionType):
        template = mypy.typeops.make_simplified_union(template.items, keep_erased=True)
    if isinstance(actual, UnionType):
        actual = mypy.typeops.make_simplified_union(actual.items, keep_erased=True)

    # Ignore Any types from the type suggestion engine to avoid them
    # causing us to infer Any in situations where a better job could
    # be done otherwise. (This can produce false positives but that
    # doesn't really matter because it is all heuristic anyway.)
    if isinstance(actual, AnyType) and actual.type_of_any == TypeOfAny.suggestion_engine:
        return []

    # If the template is simply a type variable, emit a Constraint directly.
    # We need to handle this case before handling Unions for two reasons:
    #  1. "T <: Union[U1, U2]" is not equivalent to "T <: U1 or T <: U2",
    #     because T can itself be a union (notably, Union[U1, U2] itself).
    #  2. "T :> Union[U1, U2]" is logically equivalent to "T :> U1 and
    #     T :> U2", but they are not equivalent to the constraint solver,
    #     which never introduces new Union types (it uses join() instead).
    if isinstance(template, TypeVarType):
        return [Constraint(template.id, direction, actual)]

    # Now handle the case of either template or actual being a Union.
    # For a Union to be a subtype of another type, every item of the Union
    # must be a subtype of that type, so concatenate the constraints.
    if direction == SUBTYPE_OF and isinstance(template, UnionType):
        res = []
        for t_item in template.items:
            res.extend(infer_constraints(t_item, actual, direction))
        return res
    if direction == SUPERTYPE_OF and isinstance(actual, UnionType):
        res = []
        for a_item in actual.items:
            res.extend(infer_constraints(orig_template, a_item, direction))
        return res

    # Now the potential subtype is known not to be a Union or a type
    # variable that we are solving for. In that case, for a Union to
    # be a supertype of the potential subtype, some item of the Union
    # must be a supertype of it.
    if direction == SUBTYPE_OF and isinstance(actual, UnionType):
        # If some of items is not a complete type, disregard that.
        items = simplify_away_incomplete_types(actual.items)
        # We infer constraints eagerly -- try to find constraints for a type
        # variable if possible. This seems to help with some real-world
        # use cases.
        return any_constraints(
            [infer_constraints_if_possible(template, a_item, direction)
             for a_item in items],
            eager=True)
    if direction == SUPERTYPE_OF and isinstance(template, UnionType):
        # When the template is a union, we are okay with leaving some
        # type variables indeterminate. This helps with some special
        # cases, though this isn't very principled.
        return any_constraints(
            [infer_constraints_if_possible(t_item, actual, direction)
             for t_item in template.items],
            eager=False)

    # Remaining cases are handled by ConstraintBuilderVisitor.
    return template.accept(ConstraintBuilderVisitor(actual, direction))


def infer_constraints_if_possible(template: Type, actual: Type,
                                  direction: int) -> Optional[List[Constraint]]:
    """Like infer_constraints, but return None if the input relation is
    known to be unsatisfiable, for example if template=List[T] and actual=int.
    (In this case infer_constraints would return [], just like it would for
    an automatically satisfied relation like template=List[T] and actual=object.)
    """
    if (direction == SUBTYPE_OF and
            not mypy.subtypes.is_subtype(erase_typevars(template), actual)):
        return None
    if (direction == SUPERTYPE_OF and
            not mypy.subtypes.is_subtype(actual, erase_typevars(template))):
        return None
    if (direction == SUPERTYPE_OF and isinstance(template, TypeVarType) and
            not mypy.subtypes.is_subtype(actual, erase_typevars(template.upper_bound))):
        # This is not caught by the above branch because of the erase_typevars() call,
        # that would return 'Any' for a type variable.
        return None
    return infer_constraints(template, actual, direction)


def select_trivial(options: Sequence[Optional[List[Constraint]]]) -> List[List[Constraint]]:
    """Select only those lists where each item is a constraint against Any."""
    res = []
    for option in options:
        if option is None:
            continue
        if all(isinstance(get_proper_type(c.target), AnyType) for c in option):
            res.append(option)
    return res


def merge_with_any(constraint: Constraint) -> Constraint:
    """Transform a constraint target into a union with given Any type."""
    target = constraint.target
    if is_union_with_any(target):
        # Do not produce redundant unions.
        return constraint
    # TODO: if we will support multiple sources Any, use this here instead.
    any_type = AnyType(TypeOfAny.implementation_artifact)
    return Constraint(
        constraint.type_var,
        constraint.op,
        UnionType.make_union([target, any_type], target.line, target.column),
    )


def any_constraints(options: List[Optional[List[Constraint]]], eager: bool) -> List[Constraint]:
    """Deduce what we can from a collection of constraint lists.

    It's a given that at least one of the lists must be satisfied. A
    None element in the list of options represents an unsatisfiable
    constraint and is ignored.  Ignore empty constraint lists if eager
    is true -- they are always trivially satisfiable.
    """
    if eager:
        valid_options = [option for option in options if option]
    else:
        valid_options = [option for option in options if option is not None]

    if not valid_options:
        return []

    if len(valid_options) == 1:
        return valid_options[0]

    if all(is_same_constraints(valid_options[0], c) for c in valid_options[1:]):
        # Multiple sets of constraints that are all the same. Just pick any one of them.
        return valid_options[0]

    if all(is_similar_constraints(valid_options[0], c) for c in valid_options[1:]):
        # All options have same structure. In this case we can merge-in trivial
        # options (i.e. those that only have Any) and try again.
        # TODO: More generally, if a given (variable, direction) pair appears in
        # every option, combine the bounds with meet/join always, not just for Any.
        trivial_options = select_trivial(valid_options)
        if trivial_options and len(trivial_options) < len(valid_options):
            merged_options = []
            for option in valid_options:
                if option in trivial_options:
                    continue
                if option is not None:
                    merged_option: Optional[List[Constraint]] = [
                        merge_with_any(c) for c in option
                    ]
                else:
                    merged_option = None
                merged_options.append(merged_option)
            return any_constraints([option for option in merged_options], eager)
    # Otherwise, there are either no valid options or multiple, inconsistent valid
    # options. Give up and deduce nothing.
    return []


def is_same_constraints(x: List[Constraint], y: List[Constraint]) -> bool:
    for c1 in x:
        if not any(is_same_constraint(c1, c2) for c2 in y):
            return False
    for c1 in y:
        if not any(is_same_constraint(c1, c2) for c2 in x):
            return False
    return True


def is_same_constraint(c1: Constraint, c2: Constraint) -> bool:
    # Ignore direction when comparing constraints against Any.
    skip_op_check = (
        isinstance(get_proper_type(c1.target), AnyType) and
        isinstance(get_proper_type(c2.target), AnyType)
    )
    return (c1.type_var == c2.type_var
            and (c1.op == c2.op or skip_op_check)
            and mypy.sametypes.is_same_type(c1.target, c2.target))


def is_similar_constraints(x: List[Constraint], y: List[Constraint]) -> bool:
    """Check that two lists of constraints have similar structure.

    This means that each list has same type variable plus direction pairs (i.e we
    ignore the target). Except for constraints where target is Any type, there
    we ignore direction as well.
    """
    return _is_similar_constraints(x, y) and _is_similar_constraints(y, x)


def _is_similar_constraints(x: List[Constraint], y: List[Constraint]) -> bool:
    """Check that every constraint in the first list has a similar one in the second.

    See docstring above for definition of similarity.
    """
    for c1 in x:
        has_similar = False
        for c2 in y:
            # Ignore direction when either constraint is against Any.
            skip_op_check = (
                isinstance(get_proper_type(c1.target), AnyType) or
                isinstance(get_proper_type(c2.target), AnyType)
            )
            if c1.type_var == c2.type_var and (c1.op == c2.op or skip_op_check):
                has_similar = True
                break
        if not has_similar:
            return False
    return True


def simplify_away_incomplete_types(types: Iterable[Type]) -> List[Type]:
    complete = [typ for typ in types if is_complete_type(typ)]
    if complete:
        return complete
    else:
        return list(types)


def is_complete_type(typ: Type) -> bool:
    """Is a type complete?

    A complete doesn't have uninhabited type components or (when not in strict
    optional mode) None components.
    """
    return typ.accept(CompleteTypeVisitor())


class CompleteTypeVisitor(TypeQuery[bool]):
    def __init__(self) -> None:
        super().__init__(all)

    def visit_uninhabited_type(self, t: UninhabitedType) -> bool:
        return False


class ConstraintBuilderVisitor(TypeVisitor[List[Constraint]]):
    """Visitor class for inferring type constraints."""

    # The type that is compared against a template
    # TODO: The value may be None. Is that actually correct?
    actual: ProperType

    def __init__(self, actual: ProperType, direction: int) -> None:
        # Direction must be SUBTYPE_OF or SUPERTYPE_OF.
        self.actual = actual
        self.direction = direction

    # Trivial leaf types

    def visit_unbound_type(self, template: UnboundType) -> List[Constraint]:
        return []

    def visit_any(self, template: AnyType) -> List[Constraint]:
        return []

    def visit_none_type(self, template: NoneType) -> List[Constraint]:
        return []

    def visit_uninhabited_type(self, template: UninhabitedType) -> List[Constraint]:
        return []

    def visit_erased_type(self, template: ErasedType) -> List[Constraint]:
        return []

    def visit_deleted_type(self, template: DeletedType) -> List[Constraint]:
        return []

    def visit_literal_type(self, template: LiteralType) -> List[Constraint]:
        return []

    # Errors

    def visit_partial_type(self, template: PartialType) -> List[Constraint]:
        # We can't do anything useful with a partial type here.
        assert False, "Internal error"

    # Non-trivial leaf type

    def visit_type_var(self, template: TypeVarType) -> List[Constraint]:
        assert False, ("Unexpected TypeVarType in ConstraintBuilderVisitor"
                       " (should have been handled in infer_constraints)")

    def visit_param_spec(self, template: ParamSpecType) -> List[Constraint]:
        # Can't infer ParamSpecs from component values (only via Callable[P, T]).
        return []

    # Non-leaf types

    def visit_instance(self, template: Instance) -> List[Constraint]:
        original_actual = actual = self.actual
        res: List[Constraint] = []
        if isinstance(actual, (CallableType, Overloaded)) and template.type.is_protocol:
            if template.type.protocol_members == ['__call__']:
                # Special case: a generic callback protocol
                if not any(mypy.sametypes.is_same_type(template, t)
                           for t in template.type.inferring):
                    template.type.inferring.append(template)
                    call = mypy.subtypes.find_member('__call__', template, actual,
                                                     is_operator=True)
                    assert call is not None
                    if mypy.subtypes.is_subtype(actual, erase_typevars(call)):
                        subres = infer_constraints(call, actual, self.direction)
                        res.extend(subres)
                    template.type.inferring.pop()
                    return res
        if isinstance(actual, CallableType) and actual.fallback is not None:
            actual = actual.fallback
        if isinstance(actual, Overloaded) and actual.fallback is not None:
            actual = actual.fallback
        if isinstance(actual, TypedDictType):
            actual = actual.as_anonymous().fallback
        if isinstance(actual, LiteralType):
            actual = actual.fallback
        if isinstance(actual, Instance):
            instance = actual
            erased = erase_typevars(template)
            assert isinstance(erased, Instance)  # type: ignore
            # We always try nominal inference if possible,
            # it is much faster than the structural one.
            if (self.direction == SUBTYPE_OF and
                    template.type.has_base(instance.type.fullname)):
                mapped = map_instance_to_supertype(template, instance.type)
                tvars = mapped.type.defn.type_vars
                # N.B: We use zip instead of indexing because the lengths might have
                # mismatches during daemon reprocessing.
                for tvar, mapped_arg, instance_arg in zip(tvars, mapped.args, instance.args):
                    # TODO: ParamSpecType
                    if isinstance(tvar, TypeVarType):
                        # The constraints for generic type parameters depend on variance.
                        # Include constraints from both directions if invariant.
                        if tvar.variance != CONTRAVARIANT:
                            res.extend(infer_constraints(
                                mapped_arg, instance_arg, self.direction))
                        if tvar.variance != COVARIANT:
                            res.extend(infer_constraints(
                                mapped_arg, instance_arg, neg_op(self.direction)))
                return res
            elif (self.direction == SUPERTYPE_OF and
                    instance.type.has_base(template.type.fullname)):
                mapped = map_instance_to_supertype(instance, template.type)
                tvars = template.type.defn.type_vars
                # N.B: We use zip instead of indexing because the lengths might have
                # mismatches during daemon reprocessing.
                for tvar, mapped_arg, template_arg in zip(tvars, mapped.args, template.args):
                    # TODO: ParamSpecType
                    if isinstance(tvar, TypeVarType):
                        # The constraints for generic type parameters depend on variance.
                        # Include constraints from both directions if invariant.
                        if tvar.variance != CONTRAVARIANT:
                            res.extend(infer_constraints(
                                template_arg, mapped_arg, self.direction))
                        if tvar.variance != COVARIANT:
                            res.extend(infer_constraints(
                                template_arg, mapped_arg, neg_op(self.direction)))
                return res
            if (template.type.is_protocol and self.direction == SUPERTYPE_OF and
                    # We avoid infinite recursion for structural subtypes by checking
                    # whether this type already appeared in the inference chain.
                    # This is a conservative way to break the inference cycles.
                    # It never produces any "false" constraints but gives up soon
                    # on purely structural inference cycles, see #3829.
                    # Note that we use is_protocol_implementation instead of is_subtype
                    # because some type may be considered a subtype of a protocol
                    # due to _promote, but still not implement the protocol.
                    not any(mypy.sametypes.is_same_type(template, t)
                            for t in template.type.inferring) and
                    mypy.subtypes.is_protocol_implementation(instance, erased)):
                template.type.inferring.append(template)
                res.extend(self.infer_constraints_from_protocol_members(
                    instance, template, original_actual, template))
                template.type.inferring.pop()
                return res
            elif (instance.type.is_protocol and self.direction == SUBTYPE_OF and
                  # We avoid infinite recursion for structural subtypes also here.
                  not any(mypy.sametypes.is_same_type(instance, i)
                          for i in instance.type.inferring) and
                  mypy.subtypes.is_protocol_implementation(erased, instance)):
                instance.type.inferring.append(instance)
                res.extend(self.infer_constraints_from_protocol_members(
                    instance, template, template, instance))
                instance.type.inferring.pop()
                return res
        if isinstance(actual, AnyType):
            return self.infer_against_any(template.args, actual)
        if (isinstance(actual, TupleType) and
            (is_named_instance(template, 'typing.Iterable') or
             is_named_instance(template, 'typing.Container') or
             is_named_instance(template, 'typing.Sequence') or
             is_named_instance(template, 'typing.Reversible'))
                and self.direction == SUPERTYPE_OF):
            for item in actual.items:
                cb = infer_constraints(template.args[0], item, SUPERTYPE_OF)
                res.extend(cb)
            return res
        elif isinstance(actual, TupleType) and self.direction == SUPERTYPE_OF:
            return infer_constraints(template,
                                     mypy.typeops.tuple_fallback(actual),
                                     self.direction)
        else:
            return []

    def infer_constraints_from_protocol_members(self,
                                                instance: Instance, template: Instance,
                                                subtype: Type, protocol: Instance,
                                                ) -> List[Constraint]:
        """Infer constraints for situations where either 'template' or 'instance' is a protocol.

        The 'protocol' is the one of two that is an instance of protocol type, 'subtype'
        is the type used to bind self during inference. Currently, we just infer constrains for
        every protocol member type (both ways for settable members).
        """
        res = []
        for member in protocol.type.protocol_members:
            inst = mypy.subtypes.find_member(member, instance, subtype)
            temp = mypy.subtypes.find_member(member, template, subtype)
            if inst is None or temp is None:
                return []  # See #11020
            # The above is safe since at this point we know that 'instance' is a subtype
            # of (erased) 'template', therefore it defines all protocol members
            res.extend(infer_constraints(temp, inst, self.direction))
            if (mypy.subtypes.IS_SETTABLE in
                    mypy.subtypes.get_member_flags(member, protocol.type)):
                # Settable members are invariant, add opposite constraints
                res.extend(infer_constraints(temp, inst, neg_op(self.direction)))
        return res

    def visit_callable_type(self, template: CallableType) -> List[Constraint]:
        if isinstance(self.actual, CallableType):
            res: List[Constraint] = []
            cactual = self.actual
            param_spec = template.param_spec()
            if param_spec is None:
                # FIX verify argument counts
                # FIX what if one of the functions is generic

                # We can't infer constraints from arguments if the template is Callable[..., T]
                # (with literal '...').
                if not template.is_ellipsis_args:
                    # The lengths should match, but don't crash (it will error elsewhere).
                    for t, a in zip(template.arg_types, cactual.arg_types):
                        # Negate direction due to function argument type contravariance.
                        res.extend(infer_constraints(t, a, neg_op(self.direction)))
            else:
                # TODO: Direction
                # TODO: Deal with arguments that come before param spec ones?
                res.append(Constraint(param_spec.id,
                                      SUBTYPE_OF,
                                      cactual.copy_modified(ret_type=NoneType())))

            template_ret_type, cactual_ret_type = template.ret_type, cactual.ret_type
            if template.type_guard is not None:
                template_ret_type = template.type_guard
            if cactual.type_guard is not None:
                cactual_ret_type = cactual.type_guard

            res.extend(infer_constraints(template_ret_type, cactual_ret_type,
                                         self.direction))
            return res
        elif isinstance(self.actual, AnyType):
            param_spec = template.param_spec()
            any_type = AnyType(TypeOfAny.from_another_any, source_any=self.actual)
            if param_spec is None:
                # FIX what if generic
                res = self.infer_against_any(template.arg_types, self.actual)
            else:
                res = [Constraint(param_spec.id,
                                  SUBTYPE_OF,
                                  callable_with_ellipsis(any_type, any_type, template.fallback))]
            res.extend(infer_constraints(template.ret_type, any_type, self.direction))
            return res
        elif isinstance(self.actual, Overloaded):
            return self.infer_against_overloaded(self.actual, template)
        elif isinstance(self.actual, TypeType):
            return infer_constraints(template.ret_type, self.actual.item, self.direction)
        elif isinstance(self.actual, Instance):
            # Instances with __call__ method defined are considered structural
            # subtypes of Callable with a compatible signature.
            call = mypy.subtypes.find_member('__call__', self.actual, self.actual,
                                             is_operator=True)
            if call:
                return infer_constraints(template, call, self.direction)
            else:
                return []
        else:
            return []

    def infer_against_overloaded(self, overloaded: Overloaded,
                                 template: CallableType) -> List[Constraint]:
        # Create constraints by matching an overloaded type against a template.
        # This is tricky to do in general. We cheat by only matching against
        # the first overload item that is callable compatible. This
        # seems to work somewhat well, but we should really use a more
        # reliable technique.
        item = find_matching_overload_item(overloaded, template)
        return infer_constraints(template, item, self.direction)

    def visit_tuple_type(self, template: TupleType) -> List[Constraint]:
        actual = self.actual
        if isinstance(actual, TupleType) and len(actual.items) == len(template.items):
            res: List[Constraint] = []
            for i in range(len(template.items)):
                res.extend(infer_constraints(template.items[i],
                                             actual.items[i],
                                             self.direction))
            return res
        elif isinstance(actual, AnyType):
            return self.infer_against_any(template.items, actual)
        else:
            return []

    def visit_typeddict_type(self, template: TypedDictType) -> List[Constraint]:
        actual = self.actual
        if isinstance(actual, TypedDictType):
            res: List[Constraint] = []
            # NOTE: Non-matching keys are ignored. Compatibility is checked
            #       elsewhere so this shouldn't be unsafe.
            for (item_name, template_item_type, actual_item_type) in template.zip(actual):
                res.extend(infer_constraints(template_item_type,
                                             actual_item_type,
                                             self.direction))
            return res
        elif isinstance(actual, AnyType):
            return self.infer_against_any(template.items.values(), actual)
        else:
            return []

    def visit_union_type(self, template: UnionType) -> List[Constraint]:
        assert False, ("Unexpected UnionType in ConstraintBuilderVisitor"
                       " (should have been handled in infer_constraints)")

    def visit_type_alias_type(self, template: TypeAliasType) -> List[Constraint]:
        assert False, "This should be never called, got {}".format(template)

    def infer_against_any(self, types: Iterable[Type], any_type: AnyType) -> List[Constraint]:
        res: List[Constraint] = []
        for t in types:
            # Note that we ignore variance and simply always use the
            # original direction. This is because for Any targets direction is
            # irrelevant in most cases, see e.g. is_same_constraint().
            res.extend(infer_constraints(t, any_type, self.direction))
        return res

    def visit_overloaded(self, template: Overloaded) -> List[Constraint]:
        res: List[Constraint] = []
        for t in template.items:
            res.extend(infer_constraints(t, self.actual, self.direction))
        return res

    def visit_type_type(self, template: TypeType) -> List[Constraint]:
        if isinstance(self.actual, CallableType):
            return infer_constraints(template.item, self.actual.ret_type, self.direction)
        elif isinstance(self.actual, Overloaded):
            return infer_constraints(template.item, self.actual.items[0].ret_type,
                                     self.direction)
        elif isinstance(self.actual, TypeType):
            return infer_constraints(template.item, self.actual.item, self.direction)
        elif isinstance(self.actual, AnyType):
            return infer_constraints(template.item, self.actual, self.direction)
        else:
            return []


def neg_op(op: int) -> int:
    """Map SubtypeOf to SupertypeOf and vice versa."""

    if op == SUBTYPE_OF:
        return SUPERTYPE_OF
    elif op == SUPERTYPE_OF:
        return SUBTYPE_OF
    else:
        raise ValueError('Invalid operator {}'.format(op))


def find_matching_overload_item(overloaded: Overloaded, template: CallableType) -> CallableType:
    """Disambiguate overload item against a template."""
    items = overloaded.items
    for item in items:
        # Return type may be indeterminate in the template, so ignore it when performing a
        # subtype check.
        if mypy.subtypes.is_callable_compatible(item, template,
                                                is_compat=mypy.subtypes.is_subtype,
                                                ignore_return=True):
            return item
    # Fall back to the first item if we can't find a match. This is totally arbitrary --
    # maybe we should just bail out at this point.
    return items[0]
