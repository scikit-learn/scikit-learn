'''
This module performs the return type inference, according to symbolic types,
   It then reorders function declarations according to the return type deps.
    * type_all generates a node -> type binding
'''

from pythran.analyses import LazynessAnalysis, StrictAliases, YieldPoints
from pythran.analyses import LocalNodeDeclarations, Immediates, RangeValues
from pythran.analyses import Ancestors
from pythran.analyses.aliases import ContainerOf
from pythran.config import cfg
from pythran.cxxtypes import TypeBuilder, ordered_set
from pythran.intrinsic import UserFunction, Class
from pythran.passmanager import ModuleAnalysis
from pythran.errors import PythranSyntaxError
from pythran.tables import operator_to_lambda, MODULES
from pythran.typing import List, Dict, Set, Tuple, NDArray, Union
from pythran.types.conversion import pytype_to_ctype, PYTYPE_TO_CTYPE_TABLE
from pythran.types.reorder import Reorder
from pythran.utils import attr_to_path, isnum, isextslice

from collections import defaultdict
from functools import reduce
import gast as ast
from itertools import islice
from copy import deepcopy
import numpy as np

alias_to_type = {
    MODULES['builtins']['int'] : int,
    MODULES['builtins']['bool'] : bool,
    MODULES['builtins']['float'] : float,
    MODULES['builtins']['str'] : str,
    MODULES['builtins']['complex'] : complex,
    MODULES['builtins']['dict'] : Dict,
    MODULES['builtins']['list'] : List,
    MODULES['builtins']['set'] : Set,
    MODULES['builtins']['tuple'] : Tuple,
    MODULES['builtins']['type'] : type,
    MODULES['builtins']['None'] : type(None),
    MODULES['numpy']['intc']: np.intc,
    MODULES['numpy']['intp']: np.intp,
    MODULES['numpy']['int64']: np.int64,
    MODULES['numpy']['int32']: np.int32,
    MODULES['numpy']['int16']: np.int16,
    MODULES['numpy']['int8']: np.int8,
    MODULES['numpy']['uintc']: np.uintc,
    MODULES['numpy']['uintp']: np.uintp,
    MODULES['numpy']['uint64']: np.uint64,
    MODULES['numpy']['uint32']: np.uint32,
    MODULES['numpy']['uint16']: np.uint16,
    MODULES['numpy']['uint8']: np.uint8,
    MODULES['numpy']['float32']: np.float32,
    MODULES['numpy']['float64']: np.float64,
    MODULES['numpy']['complex64']: np.complex64,
    MODULES['numpy']['complex128']: np.complex128,
    MODULES['numpy']['ndarray']: NDArray,
}
try:
    alias_to_type[MODULES['numpy']['float128']] = np.float128
    alias_to_type[MODULES['numpy']['complex256']] = np.complex256
except AttributeError:
    pass

def alias_key(a):
    if hasattr(a, 'path'):
        return a.path
    if hasattr(a, 'name'):
        return a.name,
    if hasattr(a, 'id'):
        return a.id,
    if isinstance(a, ast.Subscript):
        return ('subscript:',) + alias_key(a.value) + alias_key(a.slice)
    if isinstance(a, ast.Attribute):
        return ('attr:', a.attr) + alias_key(a.value)
    if isinstance(a, ast.Call):
        return ('call:',) + alias_key(a.func)
    if isinstance(a, ContainerOf):
        return sum((alias_key(c) for c in sorted(a.containees, key=alias_key)), (str(a.index),))
    if isinstance(a, ast.Constant):
        return ('cst:', str(a.value))
    # FIXME: how could we order those?
    return str(id(a)),



class TypeAnnotationParser(ast.NodeVisitor):

    class TypeOf:
        def __init__(self, val):
            self.val = val

    class UnionOf:
        def __init__(self, lhs, rhs):
            self.lhs = lhs
            self.rhs = rhs

    def __init__(self, type_visitor):
        self.type_visitor = type_visitor
        self.aliases = self.type_visitor.strict_aliases

    def extract(self, node):
        node_aliases = self.aliases[node]
        if len(node_aliases) > 1:
            raise PythranSyntaxError("Ambiguous identifier in type annotation",
                                     node)
        if not node_aliases:
            raise PythranSyntaxError("Unbound identifier in type annotation",
                                     node)
        node_alias, = node_aliases
        if node_alias not in alias_to_type:
            raise PythranSyntaxError("Unsupported identifier in type annotation",
                                     node)

        return alias_to_type[node_alias]


    def visit_Attribute(self, node):
        return self.extract(node)

    def visit_Name(self, node):
        return self.extract(node)

    def visit_Constant(self, node):
        return node.value

    def visit_Call(self, node):
        func = self.visit(node.func)
        if func is not type:
            raise PythranSyntaxError("Expecting a type or a call to `type(...)`",
                                     node)
        if len(node.args) != 1:
            raise PythranSyntaxError("`type` only supports a single argument",
                                     node)
        self.type_visitor.visit(node.args[0])
        return self.TypeOf(self.type_visitor.result[node.args[0]])

    def visit_Tuple(self, node):
        return tuple([self.visit(elt) for elt in node.elts])

    def visit_BinOp(self, node):
        if not isinstance(node.op, ast.BitOr):
            raise PythranSyntaxError("Unsupported operation between type operands",
                                     node)
        left = self.visit(node.left)
        right = self.visit(node.right)
        return self.UnionOf(left, right)

    def visit_Subscript(self, node):
        value = self.visit(node.value)
        slice_ = self.visit(node.slice)
        if issubclass(value, NDArray):
            dtype, ndims = slice_
            return value[tuple([dtype, *([slice(0)] * ndims)])]
        else:
            return value[slice_]


def build_type(builder, t):
    """ Python -> pythonic type binding. """
    if t is None:
        return builder.NamedType(pytype_to_ctype(type(None)))
    elif isinstance(t, List):
        return builder.ListType(build_type(builder, t.__args__[0]))
    elif isinstance(t, Set):
        return builder.SetType(build_type(builder, t.__args__[0]))
    elif isinstance(t, Dict):
        tkey, tvalue = t.__args__
        return builder.DictType(build_type(builder, tkey), build_type(builder,
                                                                      tvalue))
    elif isinstance(t, Tuple):
        return builder.TupleType(*[build_type(builder, p) for p in t.__args__])
    elif isinstance(t, NDArray):
        return builder.NDArrayType(build_type(builder, t.__args__[0]), len(t.__args__) - 1)
    elif isinstance(t, TypeAnnotationParser.TypeOf):
        return t.val
    elif isinstance(t, TypeAnnotationParser.UnionOf):
        return builder.CombinedTypes(build_type(builder, t.lhs),
                                     build_type(builder, t.rhs))
    elif t in PYTYPE_TO_CTYPE_TABLE:
        return builder.NamedType(PYTYPE_TO_CTYPE_TABLE[t])
    else:
        raise NotImplementedError("build_type on {}".format(type(t)))


def parse_type_annotation(type_visitor, ann):
    tap = TypeAnnotationParser(type_visitor)
    typ = tap.visit(ann)
    return build_type(type_visitor.builder, typ)


class UnboundableRValue(Exception):
    pass


class Types(ModuleAnalysis[Reorder, StrictAliases, LazynessAnalysis,
                           Immediates, RangeValues, Ancestors]):

    """ Infer symbolic type for all AST node. """
    class ResultType(dict):
        def __init__(self):
            self.builder = TypeBuilder()

        def copy(self):
            other = TypeResult()
            other.update(self.items())
            other.builder = self.builder
            return other


    def __init__(self):
        super().__init__()
        self.max_seq_size = cfg.getint('typing',
                                       'max_heterogeneous_sequence_size')
        self.builder = self.result.builder
        self.result["bool"] = self.builder.NamedType("bool")
        self.combiners = defaultdict(UserFunction)
        self.current_global_declarations = dict()
        self.max_recompute = 1  # max number of use to be lazy
        self.curr_locals_declaration = None
        self.ptype_count = 0

    def combined(self, *types):
        all_types = ordered_set()
        for ty in types:
            if isinstance(ty, self.builder.CombinedTypes):
                all_types.extend(ty.types)
            elif ty is not self.builder.UnknownType:
                all_types.append(ty)

        if len(all_types) == 0:
            return self.builder.UnknownType
        elif len(all_types) == 1:
            return next(iter(all_types))
        else:
            all_types = all_types[:cfg.getint('typing', 'max_combiner')]
            if {type(ty) for ty in all_types} == {self.builder.ListType}:
                return self.builder.ListType(self.combined(*[ty.of for ty in all_types]))
            if {type(ty) for ty in all_types} == {self.builder.SetType}:
                return self.builder.SetType(self.combined(*[ty.of for ty in all_types]))
            if {type(ty) for ty in all_types} == {self.builder.Assignable}:
                return self.builder.Assignable(self.combined(*[ty.of for ty in all_types]))
            if {type(ty) for ty in all_types} == {self.builder.Lazy}:
                return self.builder.Lazy(self.combined(*[ty.of for ty in all_types]))

            return self.builder.CombinedTypes(*all_types)



    def prepare(self, node):
        """
        Initialise values to prepare typing computation.

        Reorder functions to avoid dependencies issues and prepare typing
        computation setting typing values for Pythonic functions.
        """

        def register(name, module):
            """ Recursively save function typing and combiners for Pythonic."""
            for fname, function in module.items():
                if isinstance(function, dict):
                    register(name + "::" + fname, function)
                else:
                    tname = 'pythonic::{0}::functor::{1}'.format(name, fname)
                    self.result[function] = self.builder.NamedType(tname)
                    self.combiners[function] = function
                    if isinstance(function, Class):
                        register(name + "::" + fname, function.fields)

        for mname, module in MODULES.items():
            register(mname, module)
        super(Types, self).prepare(node)

    def visit_Module(self, node):
        self.generic_visit(node)
        for head in self.current_global_declarations.values():
            if head not in self.result:
                self.result[head] = "pythonic::types::none_type"

    def register(self, fname, nid, ptype):
        """register ptype as a local typedef"""
        # Too many of them leads to memory burst
        if len(self.typedefs[fname, nid]) < cfg.getint('typing',
                                                       'max_interprocedural_combiner'):
            self.typedefs[fname, nid].append(ptype)
            return True
        return False

    def node_to_id(self, n, depth=()):
        name, depth = self.node_to_name(n, depth)
        return name.id, depth

    def node_to_name(self, n, depth=()):
        if isinstance(n, ast.Name):
            return (n, depth)
        elif isinstance(n, ast.Subscript):
            if isinstance(n.slice, ast.Slice):
                return self.node_to_name(n.value, depth)
            else:
                index = n.slice.value if isnum(n.slice) else None
                return self.node_to_name(n.value, depth + (index,))
        # use alias information if any
        elif isinstance(n, ast.Call):
            for alias in self.sorted_strict_aliases(n):
                if alias is n:  # no specific alias info
                    continue
                try:
                    return self.node_to_name(alias, depth)
                except UnboundableRValue:
                    continue
        raise UnboundableRValue()

    def isargument(self, node):
        """ checks whether node aliases to a parameter."""
        try:
            node_id, _ = self.node_to_id(node)
            return (node_id in self.name_to_nodes and
                    any(isinstance(n, ast.Name) and
                        isinstance(n.ctx, ast.Param)
                        for n in self.name_to_nodes[node_id]))
        except UnboundableRValue:
            return False

    def combine(self, node, op, othernode):
        """
        Change `node` typing with combination of `node` and `othernode`.
        """
        if self.result[othernode] is self.builder.UnknownType:
            if node not in self.result:
                self.result[node] = self.builder.UnknownType
            return

        if op is None:
            op = lambda x: x

        node_aliases = ordered_set([node])
        if node in self.strict_aliases:
            node_aliases.extend(self.sorted_strict_aliases(node))
        for a in node_aliases:
            self.combine_(a, op, othernode)

    def combine_(self, node, op, othernode):
        # This comes from an assignment,so we must check where the value is
        # assigned
        name = None
        try:
            name, depth = self.node_to_name(node)
            if depth:
                former_op = op

                # update the type to reflect container nesting
                def merge_container_type(ty, index):
                    # integral index make it possible to correctly
                    # update tuple type
                    if isinstance(index, int):
                        kty = self.builder.IntegralConstant(
                                self.builder.NamedType("long"),
                                index)
                        return self.builder.IndexableContainerType(kty,
                                                                   ty)
                    elif isinstance(index, float):
                        kty = self.builder.NamedType('double')
                        return self.builder.IndexableContainerType(kty,
                                                                   ty)
                    else:
                        # FIXME: what about other key types?
                        return self.builder.ContainerType(ty)

                for node_alias in self.sorted_strict_aliases(name,
                                                             extra=[name]):
                    def traverse_alias(alias, depth, l):
                        # not interested in aliases too deep
                        if len(depth) <= l:
                            return
                        if isinstance(alias, ContainerOf):
                            d = depth[l]
                            if d is None or np.isnan(alias.index) or d == alias.index:
                                for containee in sorted(alias.containees, key=alias_key):
                                    traverse_alias(containee, depth, l + 1)
                        else:
                            def local_op(*args):
                                t = reduce(merge_container_type,
                                              depth[l:],
                                              former_op(*args))
                                return t
                            self.combine_(alias, local_op, othernode)

                    traverse_alias(node_alias, depth, 0)
                return
        except UnboundableRValue:
            pass

        if isinstance(node, ContainerOf):
            def containeeop(*args):
                container_type = op(*args)
                if isinstance(container_type, self.builder.IndexableType):
                    raise NotImplementedError
                if isinstance(container_type, (self.builder.ListType,
                                               self.builder.SetType)):
                    return container_type.of
                if isinstance(container_type, self.builder.IndexableContainerType):
                    if np.isnan(node.index):
                        return container_type.of_val
                    elif isinstance(container_type.of_key, self.builder.IntegralConstant):
                        if container_type.of_key.index != node.index:
                            raise NotImplementedError
                return self.builder.ElementType(
                        0 if np.isnan(node.index) else node.index,
                        container_type)

            for containee in sorted(node.containees, key=alias_key):
                try:
                    self.combine(containee, containeeop, othernode)
                except NotImplementedError:
                    pass

        # perform inter procedural combination
        if self.isargument(node):
            node_id, _ = self.node_to_id(node)
            if node not in self.result:
                self.result[node] = op(self.result[othernode])
            self.name_to_nodes[name.id].append(node)
            assert self.result[node], "found an alias with a type"

            parametric_type = self.builder.PType(self.current,
                                                 self.result[othernode],
                                                 self.ptype_count)
            self.ptype_count += 1

            if self.register(self.current, node_id, parametric_type):

                current_function = self.combiners[self.current]

                def translator_generator(args, op):
                    ''' capture args for translator generation'''
                    def interprocedural_type_translator(s, n):
                        translated_othernode = ast.Name(
                            '__fake__', ast.Load(), None, None)
                        s.result[translated_othernode] = (
                            parametric_type.instanciate(
                                s.current,
                                [s.result[arg] for arg in n.args]))

                        # look for modified argument
                        for p, effective_arg in enumerate(n.args):
                            formal_arg = args[p]
                            if formal_arg.id == node_id:
                                translated_node = effective_arg
                                break
                        try:
                            s.combine(translated_node,
                                      op,
                                      translated_othernode)
                        except NotImplementedError:
                            pass
                            # this may fail when the effective
                            # parameter is an expression
                        except UnboundLocalError:
                            pass
                            # this may fail when translated_node
                            # is a default parameter
                    return interprocedural_type_translator

                translator = translator_generator(self.current.args.args, op)
                current_function.add_combiner(translator)
        else:
            self.update_type(node, op, self.result[othernode])
            if name is not None:
                self.name_to_nodes[name.id].append(node)

    def update_type(self, node, ty_builder, *args):
        if ty_builder is None:
            ty, = args
        # propagate UnknowType status if one of *args is itself unknown.
        elif any(arg is self.builder.UnknownType for arg in args):
            ty = self.builder.UnknownType
        else:
            ty = ty_builder(*args)

        curr_ty = self.result.get(node, self.builder.UnknownType)
        if isinstance(curr_ty, tuple):
            return

        self.result[node] = self.combined(curr_ty, ty)

    def visit_FunctionDef(self, node):
        self.delayed_nodes = set()
        self.curr_locals_declaration = self.gather(
            LocalNodeDeclarations,
            node)
        self.current = node
        self.typedefs = defaultdict(list)
        self.name_to_nodes = defaultdict(ordered_set)
        for arg in node.args.args:
            self.name_to_nodes[arg.id].append(arg)

        self.yield_points = self.gather(YieldPoints, node)

        self.generic_visit(node)

        for delayed_node in self.delayed_nodes:
            delayed_type = self.result[delayed_node]
            if not isinstance(delayed_type, self.builder.LType):
                continue
            all_types = ordered_set(self.result[ty] for ty in
                                    self.name_to_nodes[delayed_node.id])
            final_type = self.combined(*all_types)
            delayed_type.final_type = final_type
            if final_type is delayed_type.orig:
                self.result[delayed_node] = delayed_type.orig

        # propagate type information through all aliases
        for name, nodes in self.name_to_nodes.items():
            all_types = ordered_set(self.result[ty] for ty in nodes)
            final_type = self.combined(*all_types)
            for n in nodes:
                if n not in self.delayed_nodes:
                    self.result[n] = final_type
        self.current_global_declarations[node.name] = node

        # return type may be unset if the function always raises
        return_type = self.result.get(
            node,
            self.builder.NamedType("pythonic::types::none_type"))

        self.result[node] = (
            self.builder.Returnable(return_type),
            reduce(list.__add__, self.typedefs.values(), []))
        for k in self.curr_locals_declaration:
            self.result[k] = self.get_qualifier(k)(self.result[k])

    def assignable(self, ty):
        if isinstance(ty, (self.builder.Assignable, self.builder.ListType,
                           self.builder.NamedType)):
            return ty
        else:
            return self.builder.Assignable(ty)

    def lazy(self, ty):
        if isinstance(ty, (self.builder.Lazy, self.builder.ListType,
                           self.builder.NamedType)):
            return ty
        else:
            return self.builder.Lazy(ty)

    def get_qualifier(self, node):
        lazy_res = self.lazyness_analysis[node.id]
        return (self.lazy
                if lazy_res <= self.max_recompute
                else self.assignable)

    def visit_Return(self, node):
        """ Compute return type and merges with others possible return type."""
        self.generic_visit(node)
        # No merge are done if the function is a generator.
        if not self.yield_points:
            assert node.value, "Values were added in each return statement."
            self.update_type(self.current, None, self.result[node.value])

    def visit_Yield(self, node):
        """ Compute yield type and merges it with others yield type. """
        self.generic_visit(node)
        self.update_type(self.current, None, self.result[node.value])

    def visit_Assign(self, node):
        self.visit(node.value)
        for t in node.targets:
            self.combine(t, None, node.value)
            if t in self.curr_locals_declaration:
                self.result[t] = self.get_qualifier(t)(self.result[t])
            if isinstance(t, ast.Subscript):
                self.visit_AssignedSubscript(t)

    def visit_AnnAssign(self, node):
        node_type = parse_type_annotation(self, node.annotation)
        t = node.target
        if node_type:
            self.result[t] = node_type
        if not node.value:
            self.curr_locals_declaration.remove(t)
            return
        self.visit(node.value)
        if node_type:
            # A bit odd, isn't it? :-)
            self.combine(t, None, t)
        else:
            self.combine(t, None, node.value)
        if t in self.curr_locals_declaration:
            self.result[t] = self.get_qualifier(t)(self.result[t])

        if isinstance(t, ast.Subscript):
            self.visit_AssignedSubscript(t)

    def visit_AugAssign(self, node):
        # No visit_AssignedSubscript as the container should already have been
        # populated.
        if isinstance(node.target, ast.Subscript):
            self.visit(node.target)
            self.visit(node.value)
        else:
            tmp = ast.BinOp(deepcopy(node.target), node.op, node.value)
            self.visit(tmp)
            self.combine(node.target, None, tmp)


    def visit_For(self, node):
        self.visit(node.iter)
        self.combine(node.target, self.builder.IteratorContentType, node.iter)
        for n in node.body + node.orelse:
            self.visit(n)

    def visit_BoolOp(self, node):
        """
        Merge BoolOp operand type.

        BoolOp are "and" and "or" and may return any of these results so all
        operands should have the combinable type.
        """
        # Visit subnodes
        self.generic_visit(node)
        # Merge all operands types.
        for value in node.values:
            self.update_type(node, None, self.result[value])

    def visit_BinOp(self, node):
        self.generic_visit(node)
        self.update_type(node,
                         self.builder.ExpressionType,
                         operator_to_lambda[type(node.op)],
                         self.result[node.left],
                         self.result[node.right])

    def visit_UnaryOp(self, node):
        self.generic_visit(node)
        self.update_type(node,
                         self.builder.ExpressionType,
                         operator_to_lambda[type(node.op)],
                         self.result[node.operand])

    def visit_IfExp(self, node):
        self.generic_visit(node)
        self.update_type(node, None, self.result[node.body])
        self.update_type(node, None, self.result[node.orelse])

    def visit_Compare(self, node):
        self.generic_visit(node)
        all_compare = zip(node.ops,
                          [node.left] + node.comparators[:-1],
                          node.comparators)

        for op, left, right in all_compare:
            self.update_type(node,
                             self.builder.ExpressionType,
                             operator_to_lambda[type(op)],
                             self.result[left],
                             self.result[right])

    def sorted_strict_aliases(self, func, extra=[]):
        return sorted(self.strict_aliases[func] | set(extra), key=alias_key)

    def visit_Call(self, node):
        self.generic_visit(node)

        func = node.func

        for alias in self.sorted_strict_aliases(func):
            # this comes from a bind
            if isinstance(alias, ast.Call):
                a0 = alias.args[0]
                # by construction of the bind construct
                assert len(self.strict_aliases[a0]) == 1
                bounded_function = next(iter(self.sorted_strict_aliases(a0)))
                fake_name = deepcopy(a0)
                fake_node = ast.Call(fake_name, alias.args[1:] + node.args,
                                     [])
                self.combiners[bounded_function].combiner(self, fake_node)

            # handle backward type dependencies from function calls
            else:
                self.combiners[alias].combiner(self, node)

        UnknownType = self.builder.UnknownType

        # recurring nightmare
        def last_chance():
            # maybe we can get saved if we have a hint about
            # the called function return type
            for alias in self.sorted_strict_aliases(func):
                if alias is self.current and alias in self.result:
                    # great we have a (partial) type information
                    self.result[node] = self.result[alias]
                    return
            self.result[node] = UnknownType

        if self.result[node.func] is UnknownType:
            return last_chance()

        if any(self.result[arg] is UnknownType for arg in node.args):
            return last_chance()

        # special handler for getattr: use the attr name as an enum member
        if (isinstance(func, ast.Attribute) and func.attr == 'getattr'):
            self.update_type(node,
                             self.builder.GetAttr,
                             self.result[node.args[0]],
                             node.args[1].value)
        # default behavior
        else:
            self.update_type(node,
                             self.builder.ReturnType,
                             self.result[func],
                             *[self.result[arg] for arg in node.args])

    def visit_Constant(self, node):
        """ Set the pythonic constant type. """
        ty = type(node.value)
        if ty is str and len(node.value) == 1:
            sty = 'pythonic::types::chr'
        else:
            sty = pytype_to_ctype(ty)
        if node in self.immediates:
            if sty == 'pythonic::types::chr':
                sty = self.builder.IntegralConstant(
                        self.builder.NamedType("char"),
                        f"'{node.value}'")
            else:
                sty = self.builder.IntegralConstant(
                        self.builder.NamedType(sty),
                        node.value)
        else:
            sty = self.builder.NamedType(sty)
        self.result[node] = sty

    def visit_Attribute(self, node):
        """ Compute typing for an attribute node. """
        obj, path = attr_to_path(node)
        # If no type is given, use a decltype
        if obj.isliteral():
            typename = pytype_to_ctype(obj.signature)
            self.result[node] = self.builder.NamedType(typename)
        else:
            self.result[node] = self.builder.FunctionType(*path)

    def visit_Slice(self, node):
        """
        Set slicing type using continuous information if provided.

        Also visit subnodes as they may contains relevant typing information.
        """
        self.generic_visit(node)
        nstep = node.step
        if nstep is None or (isnum(nstep) and nstep.value > 0):
            if nstep is None or nstep.value == 1:
                if all(self.range_values[p].low >= 0
                       for p in (node.lower, node.upper)):
                    ntype = "pythonic::types::fast_contiguous_slice"
                else:
                    ntype = "pythonic::types::contiguous_slice"
            else:
                ntype = "pythonic::types::cstride_slice<{}>".format(nstep.value)
            self.result[node] = self.builder.NamedType(ntype)
        else:
            self.result[node] = self.builder.NamedType(
                'pythonic::types::slice')

    def stores_to(self, node):
        ancestors = self.ancestors[node] + (node,)
        stmt_indices = [i for i, n in enumerate(ancestors)
                        if isinstance(n, (ast.Assign, ast.For))]
        if not stmt_indices:
            return True

        stmt_index = stmt_indices[-1]

        if isinstance(ancestors[stmt_index], ast.Assign):
            return ancestors[stmt_index + 1] is ancestors[stmt_index].value
        else:
            return ancestors[stmt_index + 1] is not ancestors[stmt_index].target

    def visit_Subscript(self, node):
        self.visit(node.value)
        if self.stores_to(node):
            self.result[node.value] = self.builder.AddConst(self.result[node.value])
        # type of a[1:2, 3, 4:1] is the type of: declval(a)(slice, long, slice)
        if isextslice(node.slice):
            self.visit(node.slice)
            self.update_type(node,
                             self.builder.ExpressionType,
                             lambda a, *b: "{0}({1})".format(a, ", ".join(b)),
                             self.result[node.value],
                             *[self.result[d] for d in node.slice.elts])
        elif isnum(node.slice) and node.slice.value >= 0:
            # type of a[2] is the type of an elements of a
            # this special case is to make type inference easier
            # for the back end compiler
            self.update_type(node,
                             self.builder.ElementType,
                             node.slice.value,
                             self.result[node.value])
        else:
            # type of a[i] is the return type of the matching function
            self.visit(node.slice)
            self.update_type(node,
                             self.builder.ExpressionType,
                             "{0}[{1}]".format,
                             self.result[node.value],
                             self.result[node.slice])

    def visit_AssignedSubscript(self, node):
        if isinstance(node.slice, ast.Slice):
            return
        elif isextslice(node.slice):
            return
        else:
            self.visit(node.slice)
            self.combine(node.value, self.builder.IndexableType, node.slice)

    def delayed(self, node):
        fallback_type = self.combined(*[self.result[n] for n in
                                        self.name_to_nodes[node.id]])
        self.delayed_nodes.add(node)
        return self.builder.LType(fallback_type, node)

    def visit_Name(self, node):
        if node.id in self.name_to_nodes:
            self.result[node] = self.delayed(node)
        elif node.id in self.current_global_declarations:
            newtype = self.builder.NamedType(
                self.current_global_declarations[node.id].name)
            if node not in self.result:
                self.result[node] = newtype
        else:
            self.result[node] = self.builder.UnknownType

    def visit_List(self, node):
        """ Define list type from all elements type (or empty_list type). """
        self.generic_visit(node)
        if node.elts:
            for elt in node.elts[:self.max_seq_size]:
                self.update_type(node, self.builder.ListType, self.result[elt])
        else:
            self.update_type(node,
                             self.builder.NamedType,
                             "pythonic::types::empty_list")

    def visit_Set(self, node):
        """ Define set type from all elements type (or empty_set type). """
        self.generic_visit(node)
        if node.elts:
            for elt in node.elts[:self.max_seq_size]:
                self.update_type(node, self.builder.SetType, self.result[elt])
        else:
            self.update_type(node, self.builder.NamedType,
                "pythonic::types::empty_set")

    def visit_Dict(self, node):
        """ Define set type from all elements type (or empty_dict type). """
        self.generic_visit(node)
        if node.keys:
            for key, value in islice(zip(node.keys, node.values),
                                     self.max_seq_size):
                self.update_type(node,
                                 self.builder.DictType,
                                 self.result[key],
                                 self.result[value])
        else:
            self.update_type(node,
                             self.builder.NamedType,
                             "pythonic::types::empty_dict")

    def visit_ExceptHandler(self, node):
        if node.type and node.name:
            if not isinstance(node.type, ast.Tuple):
                tname = self.builder.NamedType(
                    'pythonic::types::{0}'.format(node.type.attr))
                self.result[node.type] = tname
                self.combine(node.name, None, node.type)
        for n in node.body:
            self.visit(n)

    def visit_Tuple(self, node):
        self.generic_visit(node)
        types = [self.result[elt] for elt in node.elts]
        self.update_type(node,
                         self.builder.TupleType,
                         *types)

    def visit_arguments(self, node):
        for i, arg in enumerate(node.args):
            self.update_type(arg, self.builder.ArgumentType, i)
        for n in node.defaults:
            self.visit(n)
