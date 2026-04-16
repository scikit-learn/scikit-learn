""" This module contains all classes used to model intrinsics behavior.  """

from pythran.conversion import to_ast
from pythran.interval import UNKNOWN_RANGE, bool_values
from pythran.types.signature import extract_combiner
from pythran.typing import Any, Union, Fun, Generator

from itertools import repeat, chain
import gast as ast


class UnboundValueType(object):
    '''
    Represents a new location, bound to no identifier
    '''


UnboundValue = UnboundValueType()

class defaultlist:
    DefaultLen = 20
    def __init__(self, content=None, *, default):
        self.content = [] if content is None else content
        self.default = default

    def __iter__(self):
        return chain(self.content, repeat(self.default(), len(self) -
                                          len(self.content)))

    def __getitem__(self, i):
        return self.content[i] if i < len(self.content) else self.default()

    def __len__(self):
        return max(len(self.content), defaultlist.DefaultLen)


class UpdateEffect(object):
    pass


class ReadEffect(object):
    pass


class ReadOnceEffect(ReadEffect):
    pass


class Intrinsic(object):

    """
    Model any Method/Function.

    Its member variables are:

    - argument_effects that describes the effect of the function on its
      argument (either UpdateEffect, ReadEffect or ReadOnceEffect)
    - global_effects that describes whether the function has side effects
    - return_alias that describes the aliasing between the return value
      and the parameters. The lambda returns an ast expression, generally
      depending on the node arguments (see dict.setdefault)
    - args that describes the name and default value of each arg, using the
      same representation as ast.FunctionDef, i.e. ast.arguments
    """

    def __init__(self, **kwargs):
        self.argument_effects = kwargs.get('argument_effects',
                                           defaultlist(default=UpdateEffect))
        self.global_effects = kwargs.get('global_effects', False)
        self.return_alias = kwargs.get('return_alias',
                                       lambda x: {UnboundValue})
        self.args = ast.arguments(
            [ast.Name(n, ast.Param(), None, None)
             for n in kwargs.get('args', [])],
            [], None,
            [ast.Name(n, ast.Param(), None, None)
             for n in kwargs.get('kwonlyargs', [])],
            [], None,
            [to_ast(d) for d in kwargs.get('defaults', [])])
        self.return_range = kwargs.get("return_range",
                                       lambda call: UNKNOWN_RANGE)
        self.return_range_content = kwargs.get("return_range_content",
                                               lambda c: UNKNOWN_RANGE)
        self.requires_blas = kwargs.get("requires_blas", False)

    def isliteral(self):
        return False

    def isstatic(self):
        'static <=> value is known at compile time'
        return False;

    def isfunction(self):
        return False

    def isstaticfunction(self):
        return False

    def ismethod(self):
        return False

    def isattribute(self):
        return False

    def isconst(self):
        return not any(
            isinstance(x, UpdateEffect) for x in self.argument_effects
            ) and not self.global_effects

    def isreadonce(self, n):
        return isinstance(self.argument_effects[n], ReadOnceEffect)

    def combiner(self, s, node):
        pass


class FunctionIntr(Intrinsic):
    def __init__(self, **kwargs):
        kwargs.setdefault('combiners', ())
        super(FunctionIntr, self).__init__(**kwargs)
        self.combiners = kwargs['combiners']
        if 'signature' in kwargs:
            self.signature = kwargs['signature']
            deduced_combiner = extract_combiner(self.signature)
            if deduced_combiner is not None:
                self.combiners += deduced_combiner,
            if 'return_range' not in kwargs:
                if isinstance(self.signature, Union):
                    if all(r.__args__[-1] is bool
                           for r in self.signature.__args__):
                        self.return_range = bool_values
                elif isinstance(self.signature, Generator):
                    if self.signature.__args__[0] is bool:
                        self.return_range = bool_values
                elif isinstance(self.signature, Fun):
                    if self.signature.__args__[-1] is bool:
                        self.return_range = bool_values
        else:
            self.signature = Any
        if 'immediate_arguments' in kwargs:
            self.immediate_arguments = kwargs['immediate_arguments']
        else:
            self.immediate_arguments = []

    def isfunction(self):
        return True

    def isstaticfunction(self):
        return True

    def add_combiner(self, _combiner):
        self.combiners += (_combiner,)

    def combiner(self, s, node):
        for comb in self.combiners:
            comb(s, node)


class UserFunction(FunctionIntr):
    def __init__(self, *combiners, **kwargs):
        kwargs['combiners'] = combiners
        super(UserFunction, self).__init__(**kwargs)


class ConstFunctionIntr(FunctionIntr):
    def __init__(self, **kwargs):
        kwargs.setdefault('argument_effects',
                          defaultlist(default=ReadEffect))
        super(ConstFunctionIntr, self).__init__(**kwargs)


class ConstExceptionIntr(ConstFunctionIntr):
    def __init__(self, **kwargs):
        kwargs.setdefault('argument_effects',
                          defaultlist(default=ReadEffect))
        super(ConstExceptionIntr, self).__init__(**kwargs)


class ReadOnceFunctionIntr(ConstFunctionIntr):
    def __init__(self, **kwargs):
        super(ReadOnceFunctionIntr, self).__init__(
            argument_effects=defaultlist(default=ReadOnceEffect), **kwargs)


class MethodIntr(FunctionIntr):
    def __init__(self, *combiners, **kwargs):
        kwargs.setdefault('argument_effects',
                          defaultlist((UpdateEffect(),), default=ReadEffect))
        kwargs['combiners'] = combiners
        super(MethodIntr, self).__init__(**kwargs)

    def ismethod(self):
        return True

    def isstaticfunction(self):
        return False


class ConstMethodIntr(MethodIntr):
    def __init__(self, *combiners, **kwargs):
        kwargs.setdefault('argument_effects', defaultlist(default=ReadEffect))
        super(ConstMethodIntr, self).__init__(*combiners, **kwargs)


class ReadOnceMethodIntr(ConstMethodIntr):
    def __init__(self, **kwargs):
        super(ReadOnceMethodIntr, self).__init__(
            argument_effects=defaultlist(default=ReadOnceEffect), **kwargs)


class AttributeIntr(Intrinsic):

    """
    Internal representation for any attributes.

    Examples
    --------
    >> a.real
    """

    def __init__(self, **kwargs):
        """ Forward arguments. """
        super(AttributeIntr, self).__init__(**kwargs)
        if 'signature' in kwargs:
            self.signature = kwargs['signature']
        else:
            self.signature = Any

    def isattribute(self):
        """ Mark this intrinsic as an attribute. """
        return True


class StaticAttributeIntr(AttributeIntr):

    """
    Internal representation for any attributes whose value is known at
    compile-time, and is an int
    """

    def isstatic(self):
        return True

class ConstantIntr(Intrinsic):

    """
    Internal representation for any constant.

    Examples
    --------
    >> math.pi
    """

    def __init__(self, **kwargs):
        """ Forward arguments and remove arguments effects. """
        kwargs["argument_effects"] = ()
        super(ConstantIntr, self).__init__(**kwargs)

    def isliteral(self):
        """ Mark this intrinsic as a literal. """
        return True


class Class(Intrinsic):
    def __init__(self, d, *args, **kwargs):
        super(Class, self).__init__(*args, **kwargs)
        self.fields = d

    def __getitem__(self, key):
        return self.fields[key]

    def __iter__(self):
        return self.fields.__iter__()

    def __contains__(self, key):
        """ Forward key content to aliased module. """
        return key in self.fields


class ClassWithReadOnceConstructor(Class, ReadOnceFunctionIntr):
    def __init__(self, d, *args, **kwargs):
        super(ClassWithReadOnceConstructor, self).__init__(d, *args, **kwargs)


class ClassWithConstConstructor(Class, ConstFunctionIntr):
    def __init__(self, d, *args, **kwargs):
        super(ClassWithConstConstructor, self).__init__(d, *args, **kwargs)


class ExceptionClass(Class, ConstExceptionIntr):
    def __init__(self, d, *args, **kwargs):
        super(ExceptionClass, self).__init__(d, *args, **kwargs)


class UFunc(Class, ConstFunctionIntr):
    """ Representation of ufunc from numpy. """
