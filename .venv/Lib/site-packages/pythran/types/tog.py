# adapted from
# http://smallshire.org.uk/sufficientlysmall/2010/04/11/\
#       a-hindley-milner-type-inference-implementation-in-python/

import gast as ast
from copy import deepcopy

import numpy
from numpy import floating, integer, complexfloating

from pythran.tables import MODULES, attributes
import pythran.typing as typing
from pythran.errors import PythranSyntaxError, PythranTypeError
from pythran.utils import isnum


class InferenceError(Exception):
    "Raised if the type inference algorithm cannot infer types successfully"


symbol_of = {
    ast.And: 'and',
    ast.Or: 'or',
    ast.Add: '+',
    ast.Sub: '-',
    ast.Mult: '*',
    ast.Div: '/',
    ast.Mod: '%',
    ast.Pow: '**',
    ast.LShift: '<<',
    ast.RShift: '>>',
    ast.BitOr: '|',
    ast.BitXor: '^',
    ast.BitAnd: '&',
    ast.FloorDiv: '//',
    ast.Invert: '~',
    ast.MatMult: '@',
    ast.Not: '!',
    ast.UAdd: '+',
    ast.USub: '-',
}

NoneType_ = type(None)

# =======================================================#
# Types and type constructors


class TypeVariable(object):
    """A type variable standing for an arbitrary type.

    All type variables have a unique id, but names are only assigned lazily,
    when required.
    """

    _cached_names = {}

    def __init__(self):
        self.instance = None
        self.name = None

    def __str__(self):
        if self.instance:
            return str(self.instance)
        else:
            return 'T{}'.format(
                TypeVariable._cached_names.setdefault(
                    self,
                    len(TypeVariable._cached_names)
                )
            )


class TypeOperator(object):
    """An n-ary type constructor which builds a new type from old"""

    def __init__(self, name, types):
        self.name = name
        self.types = types

    def __str__(self):
        num_types = len(self.types)
        if num_types == 0:
            return self.name
        elif self.name == 'fun':
            return 'Callable[[{0}], {1}]'.format(
                ', '.join(map(str, self.types[:-1])), self.types[-1])
        elif self.name == 'option':
            return 'Option[{0}]'.format(self.types[0])
        else:
            return "{0}[{1}]" .format(self.name.capitalize(),
                                      ', '.join(map(str, self.types)))


class Collection(TypeOperator):

    def __init__(self, holder_type, key_type, value_type, iter_type):
        super(Collection, self).__init__("collection",
                                         [holder_type, key_type, value_type,
                                          iter_type])

    def __str__(self):
        t0 = prune(self.types[0])
        if isinstance(t0, TypeVariable):
            if isinstance(prune(self.types[1]), TypeVariable):
                return 'Iterable[{}]'.format(self.types[3])
            else:
                return 'Collection[{}, {}]'.format(self.types[1],
                                                   self.types[2])
        if isinstance(t0, TypeOperator) and t0.name == 'traits':
            if all(isinstance(prune(t), TypeVariable) for t in t0.types):
                return 'Collection[{}, {}]'.format(self.types[1],
                                                   self.types[2])
            elif all(isinstance(prune(t), TypeVariable)
                     for t in t0.types[:1] + t0.types[2:]):
                t01 = prune(t0.types[1])
                if isinstance(t01, TypeOperator) and t01.name == LenTrait.name:
                    return 'Sized'
        t00 = prune(t0.types[0])
        if isinstance(t00, TypeOperator):
            type_trait = t00.name
            if type_trait == 'list':
                return 'List[{}]'.format(self.types[2])
            if type_trait == 'set':
                return 'Set[{}]'.format(self.types[2])
            if type_trait == 'dict':
                return 'Dict[{}, {}]'.format(self.types[1], self.types[2])
            if type_trait == 'str':
                return 'str'
            if type_trait == 'file':
                return 'IO[str]'
            if type_trait == 'tuple':
                return 'Tuple[{}]'.format(', '.join(map(str, self.types[1:])))
            if type_trait == 'array':
                t01 = prune(t0.types[1])
                hasnolen = (isinstance(t01, TypeOperator) and
                            t01.name == NoLenTrait.name)
                if hasnolen:
                    return str(self.types[2])

                def rec(n):
                    pn = prune(n)
                    if isinstance(pn, Collection):
                        traits = prune(pn.types[0])
                        # a scalar or array?
                        if isinstance(traits, TypeVariable):
                            return pn.types[3], 0
                        len_trait = prune(traits.types[1])
                        # an array?
                        haslen = (isinstance(len_trait, TypeOperator) and
                                  len_trait.name == LenTrait.name)
                        if haslen:
                            t, n = rec(pn.types[3])
                            return t, n + 1
                        # a scalar or array?
                        else:
                            return pn.types[2], 0
                    else:
                        return pn, 0
                t, n = rec(self)
                if isinstance(t, TypeVariable):
                    return 'Array[{} d+, {}]'.format(n, t)
                else:
                    return 'Array[{}d, {}]'.format(n, t)
            if type_trait == 'gen':
                return 'Generator[{}]'.format(self.types[2])
        return super(Collection, self).__str__()


def TupleTrait(of_types):
    return TypeOperator('tuple', of_types)


ListTrait = TypeOperator('list', [])
SetTrait = TypeOperator('set', [])
DictTrait = TypeOperator('dict', [])
StrTrait = TypeOperator('str', [])
FileTrait = TypeOperator('file', [])
ArrayTrait = TypeOperator('array', [])
GenerableTrait = TypeOperator('gen', [])

LenTrait = TypeOperator("len", [])
NoLenTrait = TypeOperator("no_len", [])
SliceTrait = TypeOperator("slice", [])
NoSliceTrait = TypeOperator("no_slice", [])


def File():
    return Collection(Traits([FileTrait, NoLenTrait, NoSliceTrait]),
                      InvalidKey, Str(), Str())


def List(of_type):
    return Collection(Traits([ListTrait, LenTrait, SliceTrait]),
                      Integer(), of_type, of_type)


def Set(of_type):
    return Collection(Traits([SetTrait, LenTrait, NoSliceTrait]),
                      InvalidKey, of_type, of_type)


def Dict(key_type, value_type):
    return Collection(Traits([DictTrait, LenTrait, NoSliceTrait]),
                      key_type, value_type, key_type)


def Str(rec=6):
    Next = Str(rec - 1) if rec else TypeVariable()
    return Collection(Traits([StrTrait, LenTrait, SliceTrait]),
                      Integer(),
                      Next,
                      Next)


def Array(of_type, dim):
    return Collection(Traits([ArrayTrait, LenTrait, SliceTrait]),
                      AnyType,
                      AnyType,
                      Array(of_type, dim - 1) if dim > 1 else of_type)


def Iterable(of_type, dim):
    return Collection(Traits([TypeVariable(), LenTrait, SliceTrait]),
                      AnyType,
                      AnyType,
                      Iterable(of_type, dim - 1) if dim > 1 else of_type)


def Generator(of_type):
    return Collection(Traits([GenerableTrait, NoLenTrait, NoSliceTrait]),
                      InvalidKey, of_type, of_type)


def Tuple(of_types):
    return Collection(Traits([TupleTrait(of_types), LenTrait, SliceTrait]),
                      Integer(), TypeVariable(), TypeVariable())


class Scalar(TypeOperator):
    def __init__(self, types=None):
        if not isinstance(types, list):
            dtype = types
            if dtype == 'complex':
                types = [ComplexTrait,    TypeVariable(),
                         TypeVariable(), TypeVariable()]
            elif dtype == 'float':
                types = [TypeVariable(),    FloatTrait,
                         TypeVariable(), TypeVariable()]
            elif dtype == 'int':
                types = [TypeVariable(),  TypeVariable(),
                         IntegerTrait, TypeVariable()]
            elif dtype == 'bool':
                types = [TypeVariable(),  TypeVariable(),
                         TypeVariable(),      BoolTrait]
            else:
                assert dtype is None
                types = [TypeVariable(), TypeVariable(),
                         TypeVariable(), TypeVariable()]
        super(Scalar, self).__init__('scalar', types)

    def __str__(self):
        if isinstance(prune(self.types[0]), TypeOperator):
            return 'complex'
        if isinstance(prune(self.types[1]), TypeOperator):
            return 'float'
        if isinstance(prune(self.types[2]), TypeOperator):
            return 'int'
        if isinstance(prune(self.types[3]), TypeOperator):
            return 'bool'
        return 'Scalar'


def Complex():
    return Collection(Traits([ArrayTrait, NoLenTrait, NoSliceTrait]),
                      InvalidKey, Scalar('complex'), InvalidKey)


def Float():
    return Collection(Traits([ArrayTrait, NoLenTrait, NoSliceTrait]),
                      InvalidKey, Scalar('float'), InvalidKey)


def Integer():
    return Collection(Traits([ArrayTrait, NoLenTrait, NoSliceTrait]),
                      InvalidKey, Scalar('int'), InvalidKey)


def Bool():
    return Collection(Traits([ArrayTrait, NoLenTrait, NoSliceTrait]),
                      InvalidKey, Scalar('bool'), InvalidKey)


def DType():
    return Collection(Traits([ArrayTrait, NoLenTrait, NoSliceTrait]),
                      InvalidKey, Scalar(), InvalidKey)


def Function(from_types, to_type):
    """A binary type constructor which builds function types"""
    return TypeOperator('fun', list(from_types) + [to_type])


def OptionType(of_type):
    return TypeOperator("option", [of_type])


def Traits(of_types):
    return TypeOperator("traits", of_types)


ExceptionType = TypeOperator("exception", [])

# Basic types are constructed with a null type constructor
IntegerTrait = TypeOperator("int", [])  # any integer
FloatTrait = TypeOperator("float", [])  # any float
ComplexTrait = TypeOperator("complex", [])
BoolTrait = TypeOperator("bool", [])
InvalidKey = TypeOperator("invalid-key", [])  # for non-indexable collection
NoneType = TypeOperator("none", [])
AnyType = TypeOperator("any", [])
InvalidType = TypeOperator("invalid-type", [])
Slice = TypeOperator("slice", [])  # slice


def is_none(t):
    pt = prune(t)
    return isinstance(pt, TypeOperator) and pt.name == "none"


def is_option_type(t):
    pt = prune(t)
    return isinstance(pt, TypeOperator) and pt.name == "option"


def maybe_array_type(t):
    pt = prune(t)
    if isinstance(pt, TypeVariable):
        return True  # maybe an array :-/
    if isinstance(pt, TypeOperator) and pt.name == "collection":
        st = prune(pt.types[0])
        if isinstance(st, TypeOperator) and st.name == "traits":
            tt = prune(st.types[0])
            if isinstance(tt, TypeVariable):
                return True  # maybe
            return isinstance(tt, TypeOperator) and tt.name == "array"
    return False


def is_test_is_none(node):
    if not isinstance(node, ast.Compare):
        return False
    left = node.left
    comparators = node.comparators
    ops = node.ops
    if len(ops) != 1:
        return False

    op = ops[0]
    if type(op) not in (ast.Is, ast.Eq):
        return False

    comparator = comparators[0]
    if not isinstance(comparator, ast.Attribute):
        return False

    return comparator.attr == 'None' and isinstance(left, ast.Name)


def is_tuple_type(t):
    pt = prune(t)
    if isinstance(pt, TypeOperator) and pt.name == "collection":
        st = prune(pt.types[0])
        if isinstance(st, TypeOperator) and st.name == "traits":
            tt = prune(st.types[0])
            return isinstance(tt, TypeOperator) and tt.name == "tuple"
    return False


def is_getattr(node):
    if not isinstance(node, ast.Call):
        return False
    if not isinstance(node.func, ast.Attribute):
        return False
    return node.func.attr == 'getattr'


class MultiType(object):
    """A binary type constructor which builds function types"""

    def __init__(self, types):
        self.name = 'multitype'
        self.types = types

    def __str__(self):
        return '\n'.join(sorted(map(str, self.types)))


def tr(t):
    def rec_tr(t, env):
        if isinstance(t, typing.TypeVar):
            if t in env:
                return env[t]
            else:
                env[t] = TypeVariable()
                return env[t]

        elif t is typing.Any:
            return TypeVariable()

        elif isinstance(t, NoneType_):
            return NoneType

        elif t in (bool, getattr(numpy, 'bool', bool)):
            return Bool()

        elif issubclass(t, slice):
            return Slice

        elif issubclass(t, (complex, complexfloating)):
            return Complex()

        elif issubclass(t, (float, floating)):
            return Float()

        elif issubclass(t, (int, integer)):
            return Integer()

        elif issubclass(t, NoneType_):
            return NoneType

        elif t is str:
            return Str()

        elif isinstance(t, typing.Generator):
            return Generator(rec_tr(t.__args__[0], env))

        elif isinstance(t, typing.List):
            return List(rec_tr(t.__args__[0], env))

        elif isinstance(t, typing.Optional):
            return OptionType(rec_tr(t.__args__[0], env))

        elif isinstance(t, typing.Set):
            return Set(rec_tr(t.__args__[0], env))

        elif isinstance(t, typing.Dict):
            return Dict(rec_tr(t.__args__[0], env), rec_tr(t.__args__[1], env))

        elif isinstance(t, typing.Tuple):
            return Tuple([rec_tr(tp, env) for tp in t.__args__])

        elif isinstance(t, typing.NDArray):
            return Array(rec_tr(t.__args__[0], env), len(t.__args__[1:]))

        elif isinstance(t, typing.Pointer):
            return Array(rec_tr(t.__args__[0], env), 1)

        elif isinstance(t, typing.Union):
            return MultiType([rec_tr(ut, env) for ut in t.__args__])

        elif t is typing.File:
            return File()

        elif isinstance(t, typing.Iterable):
            return Collection(TypeVariable(), TypeVariable(), TypeVariable(),
                              rec_tr(t.__args__[0], env))

        elif t is typing.Sized:
            return Collection(
                Traits([TypeVariable(), LenTrait, TypeVariable()]),
                TypeVariable(), TypeVariable(), TypeVariable()
            )

        elif isinstance(t, typing.Fun):
            return Function([rec_tr(at, env) for at in t.__args__[:-1]],
                            rec_tr(t.__args__[-1], env))

        else:
            raise NotImplementedError(t)

    if isinstance(t, dict):
        return t
    elif hasattr(t, 'signature'):
        return rec_tr(t.signature, {})
    else:
        return rec_tr(t, {})


####


def analyse_body(body, env, non_generic):
    # first step to gather global symbols
    for stmt in body:
        if isinstance(stmt, ast.FunctionDef):
            new_type = TypeVariable()
            env[stmt.name] = new_type
    # second to perform local inference
    for stmt in body:
        analyse(stmt, env, non_generic)


class HasYield(ast.NodeVisitor):

    def __init__(self):
        super(HasYield, self).__init__()
        self.has_yield = False

    def visit_FunctionDef(self, node):
        pass

    def visit_Yield(self, node):
        self.has_yield = True


def analyse(node, env, non_generic=None):
    """Computes the type of the expression given by node.

    The type of the node is computed in the context of the context of the
    supplied type environment env. Data types can be introduced into the
    language simply by having a predefined set of identifiers in the initial
    environment. Environment; this way there is no need to change the syntax
    or more importantly, the type-checking program when extending the language.

    Args:
        node: The root of the abstract syntax tree.
        env: The type environment is a mapping of expression identifier names
            to type assignments.
        non_generic: A set of non-generic variables, or None

    Returns:
        The computed type of the expression.

    Raises:
        InferenceError: The type of the expression could not be inferred,
        PythranTypeError: InferenceError with user friendly message + location
    """

    if non_generic is None:
        non_generic = set()

    # expr
    if isinstance(node, ast.Name):
        if isinstance(node.ctx, (ast.Store)):
            new_type = TypeVariable()
            non_generic.add(new_type)
            env[node.id] = new_type
        return get_type(node.id, env, non_generic)
    elif isinstance(node, ast.Constant):
        if isinstance(node.value, str):
            return Str()
        elif isinstance(node.value, int):
            return Integer()
        elif isinstance(node.value, float):
            return Float()
        elif isinstance(node.value, complex):
            return Complex()
        elif node.value is None:
            return NoneType
        else:
            raise NotImplementedError
    elif isinstance(node, ast.Compare):
        left_type = analyse(node.left, env, non_generic)
        comparators_type = [analyse(comparator, env, non_generic)
                            for comparator in node.comparators]
        ops_type = [analyse(op, env, non_generic)
                    for op in node.ops]
        prev_type = left_type
        result_type = TypeVariable()
        for op_type, comparator_type in zip(ops_type, comparators_type):
            try:
                unify(Function([prev_type, comparator_type], result_type),
                      op_type)
                prev_type = comparator_type
            except InferenceError:
                raise PythranTypeError(
                    "Invalid comparison, between `{}` and `{}`".format(
                        prev_type,
                        comparator_type
                    ),
                    node)
        return result_type
    elif isinstance(node, ast.Call):
        if is_getattr(node):
            self_type = analyse(node.args[0], env, non_generic)
            attr_name = node.args[1].value
            _, attr_signature = attributes[attr_name]
            attr_type = tr(attr_signature)
            result_type = TypeVariable()
            try:
                unify(Function([self_type], result_type), attr_type)
            except InferenceError:
                if isinstance(prune(attr_type), MultiType):
                    msg = 'no attribute found, tried:\n{}'.format(attr_type)
                else:
                    msg = 'tried {}'.format(attr_type)
                raise PythranTypeError(
                    "Invalid attribute for getattr call with self"
                    "of type `{}`, {}".format(self_type, msg), node)

        else:
            fun_type = analyse(node.func, env, non_generic)
            arg_types = [analyse(arg, env, non_generic) for arg in node.args]
            result_type = TypeVariable()
            try:
                unify(Function(arg_types, result_type), fun_type)
            except InferenceError:
                # recover original type
                fun_type = analyse(node.func, env, non_generic)
                if isinstance(prune(fun_type), MultiType):
                    msg = 'no overload found, tried:\n{}'.format(fun_type)
                else:
                    msg = 'tried {}'.format(fun_type)
                raise PythranTypeError(
                    "Invalid argument type for function call to "
                    "`Callable[[{}], ...]`, {}"
                    .format(', '.join('{}'.format(at) for at in arg_types),
                            msg),
                    node)
        return result_type

    elif isinstance(node, ast.IfExp):
        test_type = analyse(node.test, env, non_generic)
        unify(Function([test_type], Bool()),
              tr(MODULES['builtins']['bool']))

        if is_test_is_none(node.test):
            none_id = node.test.left.id
            body_env = env.copy()
            body_env[none_id] = NoneType
        else:
            none_id = None
            body_env = env

        body_type = analyse(node.body, body_env, non_generic)

        if none_id:
            orelse_env = env.copy()
            if is_option_type(env[none_id]):
                orelse_env[none_id] = prune(env[none_id]).types[0]
            else:
                orelse_env[none_id] = TypeVariable()
        else:
            orelse_env = env

        orelse_type = analyse(node.orelse, orelse_env, non_generic)

        try:
            return merge_unify(body_type, orelse_type)
        except InferenceError:
            raise PythranTypeError(
                "Incompatible types from different branches:"
                "`{}` and `{}`".format(
                    body_type,
                    orelse_type
                ),
                node
            )
    elif isinstance(node, ast.UnaryOp):
        operand_type = analyse(node.operand, env, non_generic)
        op_type = analyse(node.op, env, non_generic)
        result_type = TypeVariable()
        try:
            unify(Function([operand_type], result_type), op_type)
            return result_type
        except InferenceError:
            raise PythranTypeError(
                "Invalid operand for `{}`: `{}`".format(
                    symbol_of[type(node.op)],
                    operand_type
                ),
                node
            )
    elif isinstance(node, ast.BinOp):
        left_type = analyse(node.left, env, non_generic)
        op_type = analyse(node.op, env, non_generic)
        right_type = analyse(node.right, env, non_generic)
        result_type = TypeVariable()
        try:
            unify(Function([left_type, right_type], result_type), op_type)
        except InferenceError:
            raise PythranTypeError(
                "Invalid operand for `{}`: `{}` and `{}`".format(
                    symbol_of[type(node.op)],
                    left_type,
                    right_type),
                node
            )
        return result_type
    elif isinstance(node, ast.Pow):
        return tr(MODULES['numpy']['power'])
    elif isinstance(node, ast.Sub):
        return tr(MODULES['operator']['sub'])
    elif isinstance(node, (ast.USub, ast.UAdd)):
        return tr(MODULES['operator']['pos'])
    elif isinstance(node, (ast.Eq, ast.NotEq, ast.Lt, ast.LtE, ast.Gt,
                           ast.GtE, ast.Is, ast.IsNot)):
        return tr(MODULES['operator']['eq'])
    elif isinstance(node, (ast.In, ast.NotIn)):
        contains_sig = tr(MODULES['operator']['contains'])
        contains_sig.types[:-1] = reversed(contains_sig.types[:-1])
        return contains_sig
    elif isinstance(node, ast.Add):
        return tr(MODULES['operator']['add'])
    elif isinstance(node, ast.Mult):
        return tr(MODULES['operator']['mul'])
    elif isinstance(node, ast.MatMult):
        return tr(MODULES['operator']['matmul'])
    elif isinstance(node, (ast.Div, ast.FloorDiv)):
        return tr(MODULES['operator']['floordiv'])
    elif isinstance(node, ast.Mod):
        return tr(MODULES['operator']['mod'])
    elif isinstance(node, (ast.LShift, ast.RShift)):
        return tr(MODULES['operator']['lshift'])
    elif isinstance(node, (ast.BitXor, ast.BitAnd, ast.BitOr)):
        return tr(MODULES['operator']['lshift'])
    elif isinstance(node, ast.List):
        new_type = TypeVariable()
        for elt in node.elts:
            elt_type = analyse(elt, env, non_generic)
            try:
                unify(new_type, elt_type)
            except InferenceError:
                raise PythranTypeError(
                    "Incompatible list element type `{}` and `{}`".format(
                        new_type, elt_type),
                    node
                )
        return List(new_type)
    elif isinstance(node, ast.Set):
        new_type = TypeVariable()
        for elt in node.elts:
            elt_type = analyse(elt, env, non_generic)
            try:
                unify(new_type, elt_type)
            except InferenceError:
                raise PythranTypeError(
                    "Incompatible set element type `{}` and `{}`".format(
                        new_type, elt_type),
                    node
                )
        return Set(new_type)
    elif isinstance(node, ast.Dict):
        new_key_type = TypeVariable()
        for key in node.keys:
            key_type = analyse(key, env, non_generic)
            try:
                unify(new_key_type, key_type)
            except InferenceError:
                raise PythranTypeError(
                    "Incompatible dict key type `{}` and `{}`".format(
                        new_key_type, key_type),
                    node
                )
        new_value_type = TypeVariable()
        for value in node.values:
            value_type = analyse(value, env, non_generic)
            try:
                unify(new_value_type, value_type)
            except InferenceError:
                raise PythranTypeError(
                    "Incompatible dict value type `{}` and `{}`".format(
                        new_value_type, value_type),
                    node
                )
        return Dict(new_key_type, new_value_type)
    elif isinstance(node, ast.Tuple):
        return Tuple([analyse(elt, env, non_generic) for elt in node.elts])
    elif isinstance(node, ast.Slice):
        def unify_int_or_none(t, name):
            try:
                unify(t, Integer())
            except InferenceError:
                try:
                    unify(t, NoneType)
                except InferenceError:
                    raise PythranTypeError(
                        "Invalid slice {} type `{}`, expecting int or None"
                        .format(name, t)
                    )
        if node.lower:
            lower_type = analyse(node.lower, env, non_generic)
            unify_int_or_none(lower_type, 'lower bound')
        else:
            lower_type = Integer()
        if node.upper:
            upper_type = analyse(node.upper, env, non_generic)
            unify_int_or_none(upper_type, 'upper bound')
        else:
            upper_type = Integer()
        if node.step:
            step_type = analyse(node.step, env, non_generic)
            unify_int_or_none(step_type, 'step')
        else:
            step_type = Integer()
        return Slice
    elif isinstance(node, ast.Subscript):
        new_type = TypeVariable()
        value_type = prune(analyse(node.value, env, non_generic))
        try:
            slice_type = prune(analyse(node.slice, env, non_generic))
        except PythranTypeError as e:
            raise PythranTypeError(e.msg, node)

        if isinstance(node.slice, ast.Tuple):
            nbslice = len(node.slice.elts)
            dtype = TypeVariable()
            try:
                unify(Array(dtype, nbslice), clone(value_type))
            except InferenceError:
                raise PythranTypeError(
                    "Dimension mismatch when slicing `{}`".format(value_type),
                    node)
            return TypeVariable()  # FIXME
        else:
            # handle tuples in a special way
            num = isnum(node.slice)
            if num and is_tuple_type(value_type):
                try:
                    unify(prune(prune(value_type.types[0]).types[0])
                          .types[node.slice.value],
                          new_type)
                    return new_type
                except IndexError:
                    raise PythranTypeError(
                        "Invalid tuple indexing, "
                        "out-of-bound index `{}` for type `{}`".format(
                            node.slice.value,
                            value_type),
                        node)
        try:
            unify(tr(MODULES['operator']['getitem']),
                  Function([value_type, slice_type], new_type))
        except InferenceError:
            raise PythranTypeError(
                "Invalid subscripting of `{}` by `{}`".format(
                    value_type,
                    slice_type),
                node)
        return new_type
        return new_type
    elif isinstance(node, ast.Attribute):
        from pythran.utils import attr_to_path
        obj, path = attr_to_path(node)
        if obj.signature is typing.Any:
            return TypeVariable()
        else:
            return tr(obj)

    # stmt
    elif isinstance(node, ast.Import):
        for alias in node.names:
            if alias.name not in MODULES:
                raise NotImplementedError("unknown module: %s " % alias.name)
            if alias.asname is None:
                target = alias.name
            else:
                target = alias.asname
            env[target] = tr(MODULES[alias.name])
        return env
    elif isinstance(node, ast.ImportFrom):
        if node.module not in MODULES:
            raise NotImplementedError("unknown module: %s" % node.module)
        for alias in node.names:
            if alias.name not in MODULES[node.module]:
                raise NotImplementedError(
                    "unknown function: %s in %s" % (alias.name, node.module))
            if alias.asname is None:
                target = alias.name
            else:
                target = alias.asname
            env[target] = tr(MODULES[node.module][alias.name])
        return env
    elif isinstance(node, ast.FunctionDef):
        ftypes = []
        for i in range(1 + len(node.args.defaults)):
            new_env = env.copy()
            new_non_generic = non_generic.copy()

            # reset return special variables
            new_env.pop('@ret', None)
            new_env.pop('@gen', None)

            hy = HasYield()
            for stmt in node.body:
                hy.visit(stmt)
            new_env['@gen'] = hy.has_yield

            arg_types = []
            istop = len(node.args.args) - i
            for arg in node.args.args[:istop]:
                arg_type = TypeVariable()
                new_env[arg.id] = arg_type
                new_non_generic.add(arg_type)
                arg_types.append(arg_type)
            for arg, expr in zip(node.args.args[istop:],
                                 node.args.defaults[-i:]):
                arg_type = analyse(expr, new_env, new_non_generic)
                new_env[arg.id] = arg_type

            analyse_body(node.body, new_env, new_non_generic)

            result_type = new_env.get('@ret', NoneType)

            if new_env['@gen']:
                result_type = Generator(result_type)

            ftype = Function(arg_types, result_type)
            ftypes.append(ftype)
        if len(ftypes) == 1:
            ftype = ftypes[0]
            env[node.name] = ftype
        else:
            env[node.name] = MultiType(ftypes)
        return env
    elif isinstance(node, ast.Module):
        analyse_body(node.body, env, non_generic)
        return env
    elif isinstance(node, (ast.Pass, ast.Break, ast.Continue)):
        return env
    elif isinstance(node, ast.Expr):
        analyse(node.value, env, non_generic)
        return env
    elif isinstance(node, ast.Delete):
        for target in node.targets:
            if isinstance(target, ast.Name):
                if target.id in env:
                    del env[target.id]
                else:
                    raise PythranTypeError(
                        "Invalid del: unbound identifier `{}`".format(
                            target.id),
                        node)
            else:
                analyse(target, env, non_generic)
        return env
    elif isinstance(node, ast.Print):
        if node.dest is not None:
            analyse(node.dest, env, non_generic)
        for value in node.values:
            analyse(value, env, non_generic)
        return env
    elif isinstance(node, ast.Assign):
        defn_type = analyse(node.value, env, non_generic)
        for target in node.targets:
            target_type = analyse(target, env, non_generic)
            try:
                unify(target_type, defn_type)
            except InferenceError:
                raise PythranTypeError(
                    "Invalid assignment from type `{}` to type `{}`".format(
                        target_type,
                        defn_type),
                    node)
        return env
    elif isinstance(node, ast.AnnAssign):
        defn_type = analyse(node.value, env, non_generic)
        target_type = analyse(node.target, env, non_generic)
        try:
            unify(target_type, defn_type)
        except InferenceError:
            raise PythranTypeError(
                "Invalid assignment from type `{}` to type `{}`".format(
                    target_type,
                    defn_type),
                node)
        return env
    elif isinstance(node, ast.AugAssign):
        # FIMXE: not optimal: evaluates type of node.value twice
        fake_target = deepcopy(node.target)
        fake_target.ctx = ast.Load()
        fake_op = ast.BinOp(fake_target, node.op, node.value)
        ast.copy_location(fake_op, node)
        res_type = analyse(fake_op, env, non_generic)
        target_type = analyse(node.target, env, non_generic)

        try:
            unify(target_type, res_type)
        except InferenceError:
            raise PythranTypeError(
                "Invalid update operand for `{}`: `{}` and `{}`".format(
                    symbol_of[type(node.op)],
                    res_type,
                    target_type
                ),
                node
            )
        return env
    elif isinstance(node, ast.Raise):
        return env  # TODO
    elif isinstance(node, ast.Return):
        if env['@gen']:
            return env

        if node.value is None:
            ret_type = NoneType
        else:
            ret_type = analyse(node.value, env, non_generic)
        if '@ret' in env:
            try:
                ret_type = merge_unify(env['@ret'], ret_type)
            except InferenceError:
                raise PythranTypeError(
                    "function may returns with incompatible types "
                    "`{}` and `{}`".format(env['@ret'], ret_type),
                    node
                )

        env['@ret'] = ret_type
        return env
    elif isinstance(node, ast.Yield):
        assert env['@gen']
        assert node.value is not None

        if node.value is None:
            ret_type = NoneType
        else:
            ret_type = analyse(node.value, env, non_generic)
        if '@ret' in env:
            try:
                ret_type = merge_unify(env['@ret'], ret_type)
            except InferenceError:
                raise PythranTypeError(
                    "function may yields incompatible types "
                    "`{}` and `{}`".format(env['@ret'], ret_type),
                    node
                )

        env['@ret'] = ret_type
        return env
    elif isinstance(node, ast.For):
        iter_type = analyse(node.iter, env, non_generic)
        target_type = analyse(node.target, env, non_generic)
        unify(Collection(TypeVariable(), TypeVariable(), TypeVariable(),
                         target_type),
              iter_type)
        analyse_body(node.body, env, non_generic)
        analyse_body(node.orelse, env, non_generic)
        return env
    elif isinstance(node, ast.If):
        test_type = analyse(node.test, env, non_generic)
        unify(Function([test_type], Bool()),
              tr(MODULES['builtins']['bool']))

        body_env = env.copy()
        body_non_generic = non_generic.copy()

        if is_test_is_none(node.test):
            none_id = node.test.left.id
            body_env[none_id] = NoneType
        else:
            none_id = None

        analyse_body(node.body, body_env, body_non_generic)

        orelse_env = env.copy()
        orelse_non_generic = non_generic.copy()

        if none_id:
            if is_option_type(env[none_id]):
                orelse_env[none_id] = prune(env[none_id]).types[0]
            else:
                orelse_env[none_id] = TypeVariable()
        analyse_body(node.orelse, orelse_env, orelse_non_generic)

        for var in body_env:
            if var not in env:
                if var in orelse_env:
                    try:
                        new_type = merge_unify(body_env[var], orelse_env[var])
                    except InferenceError:
                        raise PythranTypeError(
                            "Incompatible types from different branches for "
                            "`{}`: `{}` and `{}`".format(
                                var,
                                body_env[var],
                                orelse_env[var]
                            ),
                            node
                        )
                else:
                    new_type = body_env[var]
                env[var] = new_type

        for var in orelse_env:
            if var not in env:
                # may not be unified by the prev loop if a del occured
                if var in body_env:
                    new_type = merge_unify(orelse_env[var], body_env[var])
                else:
                    new_type = orelse_env[var]
                env[var] = new_type

        if none_id:
            try:
                new_type = merge_unify(body_env[none_id], orelse_env[none_id])
            except InferenceError:
                msg = ("Inconsistent types while merging values of `{}` from "
                       "conditional branches: `{}` and `{}`")
                err = msg.format(none_id,
                                 body_env[none_id],
                                 orelse_env[none_id])
                raise PythranTypeError(err, node)
            env[none_id] = new_type

        return env
    elif isinstance(node, ast.While):
        test_type = analyse(node.test, env, non_generic)
        unify(Function([test_type], Bool()),
              tr(MODULES['builtins']['bool']))

        analyse_body(node.body, env, non_generic)
        analyse_body(node.orelse, env, non_generic)
        return env
    elif isinstance(node, ast.Try):
        analyse_body(node.body, env, non_generic)
        for handler in node.handlers:
            analyse(handler, env, non_generic)
        analyse_body(node.orelse, env, non_generic)
        analyse_body(node.finalbody, env, non_generic)
        return env
    elif isinstance(node, ast.ExceptHandler):
        if(node.name):
            new_type = ExceptionType
            non_generic.add(new_type)
            if node.name.id in env:
                unify(env[node.name.id], new_type)
            else:
                env[node.name.id] = new_type
        analyse_body(node.body, env, non_generic)
        return env
    elif isinstance(node, ast.Assert):
        if node.msg:
            analyse(node.msg, env, non_generic)
        analyse(node.test, env, non_generic)
        return env
    elif isinstance(node, ast.UnaryOp):
        operand_type = analyse(node.operand, env, non_generic)
        return_type = TypeVariable()
        op_type = analyse(node.op, env, non_generic)
        unify(Function([operand_type], return_type), op_type)
        return return_type
    elif isinstance(node, ast.Invert):
        return MultiType([Function([Bool()], Integer()),
                          Function([Integer()], Integer())])
    elif isinstance(node, ast.Not):
        return tr(MODULES['builtins']['bool'])
    elif isinstance(node, ast.BoolOp):
        op_type = analyse(node.op, env, non_generic)
        value_types = [analyse(value, env, non_generic)
                       for value in node.values]

        for value_type in value_types:
            unify(Function([value_type], Bool()),
                  tr(MODULES['builtins']['bool']))

        return_type = TypeVariable()
        prev_type = value_types[0]
        for value_type in value_types[1:]:
            unify(Function([prev_type, value_type], return_type), op_type)
            prev_type = value_type
        return return_type
    elif isinstance(node, (ast.And, ast.Or)):
        x_type = TypeVariable()
        return MultiType([
            Function([x_type, x_type], x_type),
            Function([TypeVariable(), TypeVariable()], TypeVariable()),
        ])

    raise RuntimeError("Unhandled syntax node {0}".format(type(node)))


def get_type(name, env, non_generic):
    """Get the type of identifier name from the type environment env.

    Args:
        name: The identifier name
        env: The type environment mapping from identifier names to types
        non_generic: A set of non-generic TypeVariables

    Raises:
        ParseError: Raised if name is an undefined symbol in the type
            environment.
    """
    if name in env:
        if isinstance(env[name], MultiType):
            return clone(env[name])
        return fresh(env[name], non_generic)
    else:
        print("W: Undefined symbol {0}".format(name))
        return TypeVariable()


def fresh(t, non_generic):
    """Makes a copy of a type expression.

    The type t is copied. The generic variables are duplicated and the
    non_generic variables are shared.

    Args:
        t: A type to be copied.
        non_generic: A set of non-generic TypeVariables
    """

    mappings = {}  # A mapping of TypeVariables to TypeVariables

    def freshrec(tp):
        p = prune(tp)
        if isinstance(p, TypeVariable):
            if is_generic(p, non_generic):
                if p not in mappings:
                    mappings[p] = TypeVariable()
                return mappings[p]
            else:
                return p
        elif isinstance(p, dict):
            return p  # module
        elif isinstance(p, Collection):
            return Collection(*[freshrec(x) for x in p.types])
        elif isinstance(p, Scalar):
            return Scalar([freshrec(x) for x in p.types])
        elif isinstance(p, TypeOperator):
            return TypeOperator(p.name, [freshrec(x) for x in p.types])
        elif isinstance(p, MultiType):
            return MultiType([freshrec(x) for x in p.types])
        else:
            assert False, "missing freshrec case {}".format(type(p))

    return freshrec(t)


def clone(t):
    if isinstance(t, MultiType):
        return MultiType([clone(tp) for tp in t.types])
    else:
        return fresh(t, {})


def unify(t1, t2):
    """Unify the two types t1 and t2.

    Makes the types t1 and t2 the same.

    Args:
        t1: The first type to be made equivalent
        t2: The second type to be be equivalent

    Returns:
        None

    Raises:
        InferenceError: Raised if the types cannot be unified.
    """

    a = prune(t1)
    b = prune(t2)
    if isinstance(a, TypeVariable):
        if a != b:
            if occurs_in_type(a, b):
                raise InferenceError("recursive unification")
            a.instance = b
    elif isinstance(b, TypeVariable):
        unify(b, a)
    elif isinstance(a, TypeOperator) and a.name == 'any':
        return
    elif isinstance(b, TypeOperator) and b.name == 'any':
        return
    elif isinstance(a, TypeOperator) and isinstance(b, TypeOperator):
        if len(a.types) != len(b.types):
            raise InferenceError("Type length differ")
        else:
            if a.name != b.name:
                raise InferenceError("Type name differ")
        try:
            for p, q in zip(a.types, b.types):
                unify(p, q)
        except InferenceError:
            raise
    elif isinstance(a, MultiType) and isinstance(b, MultiType):
        if len(a.types) != len(b.types):
            raise InferenceError("Type lenght differ")
        for p, q in zip(a.types, b.types):
            unify(p, q)
    elif isinstance(b, MultiType):
        return unify(b, a)
    elif isinstance(a, MultiType):
        types = []
        for t in a.types:
            try:
                t_clone = fresh(t, {})
                b_clone = fresh(b, {})
                unify(t_clone, b_clone)
                types.append(t)
            except InferenceError:
                pass
        if types:
            if len(types) == 1:
                unify(clone(types[0]), b)
            else:
                # too many overloads are found,
                # so extract as many information as we can,
                # and leave the remaining over-approximated
                def try_unify(t, ts):
                    if isinstance(t, TypeVariable):
                        return
                    if any(isinstance(tp, TypeVariable) for tp in ts):
                        return
                    if any(len(tp.types) != len(t.types) for tp in ts):
                        return
                    for i, tt in enumerate(t.types):
                        its = [prune(tp.types[i]) for tp in ts]
                        if any(isinstance(it, TypeVariable) for it in its):
                            continue
                        it0 = its[0]
                        it0ntypes = len(it0.types)
                        if all(((it.name == it0.name) and
                                (len(it.types) == it0ntypes))
                               for it in its):
                            ntypes = [TypeVariable() for _ in range(it0ntypes)]
                            new_tt = TypeOperator(it0.name, ntypes)
                            new_tt.__class__ = it0.__class__
                            unify(tt, new_tt)
                            try_unify(prune(tt), [prune(it) for it in its])
                try_unify(b, types)
        else:
            raise InferenceError("No overload")
    else:
        raise RuntimeError("Not unified {} and {}".format(type(a), type(b)))


def merge_unify(t1, t2):
    p1 = prune(t1)
    p2 = prune(t2)
    if is_none(p1) and is_none(p2):
        return p1
    if is_none(p1):
        if is_option_type(p2):
            return p2
        else:
            return OptionType(p2)
    if is_none(p2):
        return merge_unify(p2, p1)
    if is_option_type(p1) and is_option_type(p2):
        unify(p1.types[0], p2.types[0])
        return p1
    if is_option_type(p1):
        unify(p1.types[0], p2)
        return p1
    if is_option_type(p2):
        return merge_unify(p2, p1)
    unify(p1, p2)
    return p1


def prune(t):
    """Returns the currently defining instance of t.

    As a side effect, collapses the list of type instances. The function Prune
    is used whenever a type expression has to be inspected: it will always
    return a type expression which is either an uninstantiated type variable or
    a type operator; i.e. it will skip instantiated variables, and will
    actually prune them from expressions to remove long chains of instantiated
    variables.

    Args:
        t: The type to be pruned

    Returns:
        An uninstantiated TypeVariable or a TypeOperator
    """
    if isinstance(t, TypeVariable):
        if t.instance is not None:
            t.instance = prune(t.instance)
            return t.instance
    return t


def is_generic(v, non_generic):
    """Checks whether a given variable occurs in a list of non-generic variables

    Note that a variables in such a list may be instantiated to a type term,
    in which case the variables contained in the type term are considered
    non-generic.

    Note: Must be called with v pre-pruned

    Args:
        v: The TypeVariable to be tested for genericity
        non_generic: A set of non-generic TypeVariables

    Returns:
        True if v is a generic variable, otherwise False
    """
    return not occurs_in(v, non_generic)


def occurs_in_type(v, type2):
    """Checks whether a type variable occurs in a type expression.

    Note: Must be called with v pre-pruned

    Args:
        v:  The TypeVariable to be tested for
        type2: The type in which to search

    Returns:
        True if v occurs in type2, otherwise False
    """
    pruned_type2 = prune(type2)
    if pruned_type2 == v:
        return True
    elif isinstance(pruned_type2, TypeOperator):
        return occurs_in(v, pruned_type2.types)
    return False


def occurs_in(t, types):
    """Checks whether a types variable occurs in any other types.

    Args:
        t:  The TypeVariable to be tested for
        types: The sequence of types in which to search

    Returns:
        True if t occurs in any of types, otherwise False
    """
    return any(occurs_in_type(t, t2) for t2 in types)


def typecheck(node):
    types = analyse(node, {'builtins': MODULES['builtins']})
    return types
