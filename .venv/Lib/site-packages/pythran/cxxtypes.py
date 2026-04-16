'''
This module defines classes needed to manipulate c++ types from pythran.
'''

from inspect import isclass
from pythran.utils import cxxid


class ordered_set(object):

    def __init__(self, elements=None):
        self.values = list()
        self.unique_values = set()
        if elements is not None:
            self.extend(elements)

    def append(self, value):
        if value not in self.unique_values:
            self.values.append(value)
            self.unique_values.add(value)

    def extend(self, values):
        for value in values:
            self.append(value)

    def __iter__(self):
        return iter(self.values)

    def __len__(self):
        return len(self.values)

    def __getitem__(self, index):
        return self.values[index]


class TypeBuilder(object):
    '''
    >>> builder = TypeBuilder()

    >>> builder.NamedType('long long')
    long long

    >>> l_ty = builder.NamedType('long')
    >>> i_ty = builder.NamedType('int')
    >>> f_ty = builder.NamedType('float')

    >>> builder.CombinedTypes(builder.NamedType('long'), builder.NamedType('char'))
    typename __combined<long,char>::type

    >>> builder.ArgumentType(4)
    std::remove_cv_t<std::remove_reference_t<argument_type4>>

    >>> builder.Assignable(builder.NamedType("long"))
    typename pythonic::assignable<long>::type

    >>> builder.Returnable(builder.NamedType("long"))
    typename pythonic::returnable<long>::type

    >>> builder.Lazy(builder.NamedType("long"))
    typename pythonic::lazy<long>::type

    >>> builder.FunctionType("toto")
    toto

    >>> builder.IteratorOfType(builder.NamedType('some'))
    typename some::iterator
    >>> builder.IteratorOfType(builder.NamedType('typename some::stuff'))
    typename some::stuff::iterator

    >>> builder.IteratorContentType(builder.NamedType('str'))
    std::remove_cv_t<typename std::iterator_traits<\
typename std::remove_reference_t<str>::iterator>::value_type>

    >>> builder.GetAttr(builder.NamedType('complex'), 'real')
    decltype(pythonic::builtins::getattr(\
pythonic::types::attr::REAL{}, std::declval<complex>()))

    >>> builder.ReturnType(builder.NamedType('math::cos'), f_ty)
    std::result_of_t<math::cos(float)>

    >>> t = builder.TupleType(i_ty, builder.NamedType('str'))
    >>> builder.ElementType(1, t)
    std::tuple_element_t<1, std::remove_reference_t<pythonic::types::make_tuple_t<int, str>>>


    >>> builder.ListType(builder.NamedType('int'))
    pythonic::types::list<std::remove_reference_t<int>>

    >>> builder.NDArrayType(builder.NamedType('int'), 1)
    pythonic::types::ndarray<int, pythonic::types::pshape<long>>

    >>> builder.SetType(builder.NamedType('int'))
    pythonic::types::set<int>

    >>> builder.TupleType(i_ty, builder.NamedType('bool'))
    pythonic::types::make_tuple_t<int, bool>

    >>> builder.DictType(builder.NamedType('int'), builder.NamedType('float'))
    pythonic::types::dict<int,float>

    >>> builder.ContainerType(builder.NamedType('int'))
    container<std::remove_reference_t<int>>

    >>> builder.IndexableType(builder.NamedType('int'))
    indexable<int>

    >>> op = lambda x,y: x + '+' + y
    >>> builder.ExpressionType(op, l_ty, i_ty)
    decltype(std::declval<long>()+std::declval<int>())
    '''

    def __init__(builder):
        builder._instances = dict()

        class Type(object):
            """
            A generic type object to be sub-classed

            The keyword arguments are used to build the internal representation
            one attribute per key with the associated value
            """

            def __new__(cls, *args, **kwargs):
                # no memoization for PType
                if cls.__name__ == 'PType':
                    return super(Type, cls).__new__(cls)

                key = cls,
                for v in args + tuple(v for k, v in sorted(kwargs.items())):
                    if isinstance(v, list):
                        v = tuple(v)
                    key += v,

                if key not in builder._instances:
                    builder._instances[key] = super(Type, cls).__new__(cls)

                return builder._instances[key]

            def __init__(self, **kwargs):
                for k, v in kwargs.items():
                    if isinstance(v, list):
                        v = tuple(v)
                    setattr(self, k, v)

            def iscombined(self):
                return False

            def sgenerate(self):
                return self.generate(lambda obj: obj.sgenerate())

            def __repr__(self):
                return self.sgenerate()  # for testing only

            def __str__(self):
                raise NotImplementedError("Call self.sgenerate() instead")

        class NamedType(Type):
            """
            A generic type object, to hold scalar types and such
            """

            def __init__(self, srepr):
                super(NamedType, self).__init__(srepr=srepr)

            def generate(self, _):
                return self.srepr

        class PType(Type):
            """
            A generic parametric type
            """

            def __init__(self, fun, ptype, index):
                super(PType, self).__init__(fun=fun,
                                            type=ptype,
                                            name=f"__ptype{index}")

            def generate(self, ctx):
                return ctx(self.type)

            def instanciate(self, caller, arguments):
                if self.fun is caller:
                    return builder.UnknownType
                else:
                    return InstantiatedType(self.fun, self.name, arguments)

        class LType(Type):
            def __init__(self, base, node):
                super(LType, self).__init__(node=node)
                self.isrec = False
                self.orig = base
                self.final_type = base

            def generate(self, ctx):
                if self.isrec:
                    return ctx(self.orig)
                else:
                    self.isrec = True
                    res = ctx(self.final_type)
                    self.isrec = False
                    return res

        class InstantiatedType(Type):
            """
            A type instantiated from a parametric type
            """
            def __init__(self, fun, name, arguments):
                super(InstantiatedType, self).__init__(fun=fun,
                                                       name=name,
                                                       arguments=arguments)

            def generate(self, ctx):
                if self.arguments:
                    args = ", ".join(ctx(arg) for arg in self.arguments)
                    template_params = f"<{args}>"
                else:
                    template_params = ""

                return f"typename {self.fun.name}::type{template_params}::{self.name}"

        class CombinedTypes(Type):
            """
            type resulting from the combination of other types

            """

            def __init__(self, *types):
                super(CombinedTypes, self).__init__(types=types)

            def iscombined(self):
                return True

            def generate(self, ctx):
                # Degenerated trees may lead to very deep call stacks.
                # In that case try hard to recover by cutting the tree
                import sys
                current_recursion_limit = sys.getrecursionlimit()
                stypes = ordered_set()
                for t in self.types:
                    try:
                        stypes.append(ctx(t))
                    except RecursionError:
                        sys.setrecursionlimit(current_recursion_limit * 2)
                        break

                if not stypes:
                    sys.setrecursionlimit(current_recursion_limit)
                    raise RecursionError
                elif len(stypes) == 1:
                    sys.setrecursionlimit(current_recursion_limit)
                    return stypes[0]
                else:
                    stmp = f'typename __combined<{",".join(stypes)}>::type'
                    sys.setrecursionlimit(current_recursion_limit)
                    return stmp

        class IntegralConstant(Type):
            """
            A generic type object, to hold scalar types and such
            """

            def __init__(self, of, index):
                super(IntegralConstant, self).__init__(of=of, index=index)

            def generate(self, ctx):
                ty = ctx(self.of)
                return f"std::integral_constant<{ty}, {str(self.index).lower()}>"

        class ArgumentType(Type):
            """
            A type to hold function arguments
            """
            def __init__(self, num):
                super(ArgumentType, self).__init__(num=num)

            def generate(self, _):
                return f'std::remove_cv_t<std::remove_reference_t<argument_type{self.num}>>'

        class DependentType(Type):
            """
            A class to be sub-classed by any type that depends on another type
            """
            def __init__(self, of):
                assert of is not None
                super(DependentType, self).__init__(of=of)

            def iscombined(self):
                return self.of.iscombined()


        class Assignable(DependentType):
            """
            A type which can be assigned

            It is used to make the difference between
            * transient types (e.g. generated from expression template)
            * assignable types (typically type of a variable)
            """

            def generate(self, ctx):
                return f'typename pythonic::assignable<{ctx(self.of)}>::type'

        class AssignableNoEscape(DependentType):
            """
            Similar to Assignable, but it doesn't escape it's declaration scope
            """

            def generate(self, ctx):
                return f'typename pythonic::assignable_noescape<{ctx(self.of)}>::type'

        class Returnable(DependentType):
            """
            A type which can be returned

            It is used to make the difference between
            * returned types (that cannot hold a reference to avoid dangling
                              reference)
            * assignable types (local to a function)

            """

            def generate(self, ctx):
                return f'typename pythonic::returnable<{ctx(self.of)}>::type'

        class Lazy(DependentType):
            """
            A type which can be a reference

            It is used to make a lazy evaluation of numpy expressions

            """

            def generate(self, ctx):
                return f'typename pythonic::lazy<{ctx(self.of)}>::type'

        class FunctionType(Type):
            """
            Gather the type of a function reference
            """

            def __init__(self, *path):
                super(FunctionType, self).__init__(path=path)

            def generate(self, _):
                return '::'.join(map(cxxid, self.path))


        class AddConst(DependentType):
            '''
            Type of an Iterator of a container
            '''
            def generate(self, ctx):
                of_type = ctx(self.of)
                return f'decltype(pythonic::types::as_const(std::declval<{of_type}>()))'

        class IteratorOfType(DependentType):
            '''
            Type of an Iterator of a container
            '''
            def generate(self, ctx):
                container_type = ctx(self.of)
                if container_type.startswith('typename'):
                    return f'{container_type}::iterator'
                else:
                    return f'typename {container_type}::iterator'

        class IteratorContentType(DependentType):
            '''
            Type of an iterator over the content of a container
            '''

            def generate(self, ctx):
                iterator_value_type = ctx(self.of)
                return f'std::remove_cv_t<typename std::iterator_traits<typename std::remove_reference_t<{iterator_value_type}>::iterator>::value_type>'

        class GetAttr(Type):
            '''
            Type of a named attribute
            '''
            def __init__(self, param, attr):
                super(GetAttr, self).__init__(param=param, attr=attr)

            def generate(self, ctx):
                attr = f'pythonic::types::attr::{self.attr.upper()}'
                return f'decltype(pythonic::builtins::getattr({attr}{{}}, std::declval<{ctx(self.param)}>()))'

        class ReturnType(Type):
            '''
            Return type of a call with arguments
            '''
            def __init__(self, ftype, *args):
                super(ReturnType, self).__init__(ftype=ftype, args=args)

            def generate(self, ctx):
                # the return type of a constructor is obvious
                cg = ctx(self.ftype)
                args = [ctx(arg) for arg in self.args]
                return f'std::result_of_t<{cg}({", ".join(args)})>'

        class ElementType(Type):
            '''
            Type of the ith element of a tuple or container
            '''

            def __init__(self, index, of):
                super(ElementType, self).__init__(of=of, index=index)

            def iscombined(self):
                return self.of.iscombined()

            def generate(self, ctx):
                return f'std::tuple_element_t<{self.index}, std::remove_reference_t<{ctx(self.of)}>>'

        class TypeType(DependentType):
            '''
            Type holding a type
            '''

            def generate(self, ctx):
                return f'pythonic::types::type_t<{ctx(self.of)}>'


        class ListType(DependentType):
            '''
            Type holding a list of stuff of the same type
            '''

            def generate(self, ctx):
                return f'pythonic::types::list<std::remove_reference_t<{ctx(self.of)}>>'

        class SetType(DependentType):
            '''
            Type holding a set of stuff of the same type
            '''

            def generate(self, ctx):
                return f'pythonic::types::set<{ctx(self.of)}>'

        class TupleType(Type):
            '''
            Type holding a tuple of stuffs of various types
            '''
            def __init__(self, *ofs):
                super(TupleType, self).__init__(ofs=ofs)

            def iscombined(self):
                return any(of.iscombined() for of in self.ofs)

            def generate(self, ctx):
                elts = (ctx(of) for of in self.ofs)
                return f'pythonic::types::make_tuple_t<{", ".join(elts)}>'

        class DictType(Type):
            '''
            Type holding a dict of stuff of the same key and value type
            '''

            def __init__(self, of_key, of_val):
                super(DictType, self).__init__(of_key=of_key, of_val=of_val)

            def iscombined(self):
                return any((of.iscombined()
                            for of in (self.of_key, self.of_val)))

            def generate(self, ctx):
                return f'pythonic::types::dict<{ctx(self.of_key)},{ctx(self.of_val)}>'

        class NDArrayType(DependentType):
            '''
            Type holding a numpy array without view
            '''
            def __init__(self, dtype, nbdims):
                super(DependentType, self).__init__(of=dtype, nbdims=nbdims)

            def generate(self, ctx):
                return f'pythonic::types::ndarray<{ctx(self.of)}, pythonic::types::pshape<{", ".join(["long"] * self.nbdims)}>>'


        class ContainerType(DependentType):
            '''
            Type of any container of stuff of the same type
            '''

            def generate(self, ctx):
                return f'container<std::remove_reference_t<{ctx(self.of)}>>'

        class IndexableType(DependentType):
            '''
            Type of any container indexed by the same type
            '''

            def generate(self, ctx):
                return f'indexable<{ctx(self.of)}>'

        class IndexableContainerType(Type):
            '''
            Type of any container of stuff of the same type,
            indexable by another type
            '''
            def __init__(self, of_key, of_val):
                super(IndexableContainerType, self).__init__(of_key=of_key,
                                                             of_val=of_val)

            def iscombined(self):
                return any((of.iscombined()
                            for of in (self.of_key, self.of_val)))

            def generate(self, ctx):
                return f'indexable_container<{ctx(self.of_key)}, std::remove_reference_t<{ctx(self.of_val)}>>'

        class ExpressionType(Type):

            """
            Result type of an operator call.
            """

            def __init__(self, op, *exprs):
                super(ExpressionType, self).__init__(op=op, exprs=exprs)

            def iscombined(self):
                return any(expr.iscombined() for expr in self.exprs)

            def generate(self, ctx):
                gexprs = [f"std::declval<{ctx(expr)}>()" for expr in self.exprs]
                return f'decltype({self.op(*gexprs)})'

        builder.UnknownType = Type()

        for objname, obj in locals().items():
            if isclass(obj):
                setattr(builder, objname, obj)
