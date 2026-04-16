import cython
cython.declare(PyrexTypes=object, Naming=object, ExprNodes=object, Nodes=object,
               Options=object, UtilNodes=object, LetNode=object,
               LetRefNode=object, TreeFragment=object, EncodedString=object,
               error=object, warning=object, copy=object, hashlib=object, sys=object,
               itemgetter=object)

import copy
import hashlib
import sys
from operator import itemgetter

from . import Code
from . import PyrexTypes
from . import Naming
from . import ExprNodes
from . import Nodes
from . import Options
from . import Builtin
from . import Errors

from .Visitor import VisitorTransform, TreeVisitor
from .Visitor import CythonTransform, EnvTransform, ScopeTrackingTransform
from .UtilNodes import LetNode, LetRefNode
from .TreeFragment import TreeFragment
from .StringEncoding import EncodedString
from .Errors import error, warning, CompileError, InternalError


class SkipDeclarations:
    """
    Variable and function declarations can often have a deep tree structure,
    and yet most transformations don't need to descend to this depth.

    Declaration nodes are removed after AnalyseDeclarationsTransform, so there
    is no need to use this for transformations after that point.
    """
    def visit_CTypeDefNode(self, node):
        return node

    def visit_CVarDefNode(self, node):
        return node

    def visit_CDeclaratorNode(self, node):
        return node

    def visit_CBaseTypeNode(self, node):
        return node

    def visit_CEnumDefNode(self, node):
        return node

    def visit_CStructOrUnionDefNode(self, node):
        return node

    def visit_CppClassNode(self, node):
        if node.visibility != "extern":
            # Need to traverse methods.
            self.visitchildren(node)
        return node


class NormalizeTree(CythonTransform):
    """
    This transform fixes up a few things after parsing
    in order to make the parse tree more suitable for
    transforms.

    a) After parsing, blocks with only one statement will
    be represented by that statement, not by a StatListNode.
    When doing transforms this is annoying and inconsistent,
    as one cannot in general remove a statement in a consistent
    way and so on. This transform wraps any single statements
    in a StatListNode containing a single statement.

    b) The PassStatNode is a noop and serves no purpose beyond
    plugging such one-statement blocks; i.e., once parsed a
`    "pass" can just as well be represented using an empty
    StatListNode. This means less special cases to worry about
    in subsequent transforms (one always checks to see if a
    StatListNode has no children to see if the block is empty).
    """

    def __init__(self, context):
        super().__init__(context)
        self.is_in_statlist = False
        self.is_in_expr = False

    def visit_ModuleNode(self, node):
        self.visitchildren(node)
        if not isinstance(node.body, Nodes.StatListNode):
            # This can happen when the body only consists of a single (unused) declaration and no statements.
            node.body = Nodes.StatListNode(pos=node.pos, stats=[node.body])
        return node

    def visit_ExprNode(self, node):
        stacktmp = self.is_in_expr
        self.is_in_expr = True
        self.visitchildren(node)
        self.is_in_expr = stacktmp
        return node

    def visit_StatNode(self, node, is_listcontainer=False):
        stacktmp = self.is_in_statlist
        self.is_in_statlist = is_listcontainer
        self.visitchildren(node)
        self.is_in_statlist = stacktmp
        if not self.is_in_statlist and not self.is_in_expr:
            return Nodes.StatListNode(pos=node.pos, stats=[node])
        else:
            return node

    def visit_StatListNode(self, node):
        self.is_in_statlist = True
        self.visitchildren(node)
        self.is_in_statlist = False
        return node

    def visit_ParallelAssignmentNode(self, node):
        return self.visit_StatNode(node, True)

    def visit_CEnumDefNode(self, node):
        return self.visit_StatNode(node, True)

    def visit_CStructOrUnionDefNode(self, node):
        return self.visit_StatNode(node, True)

    def visit_ExprStatNode(self, node):
        """Eliminate useless string literals"""
        if node.expr.is_string_literal:
            return Nodes.PassStatNode(node.expr.pos)
        else:
            return self.visit_StatNode(node)

    def visit_CDeclaratorNode(self, node):
        return node


class PostParseError(CompileError): pass

# error strings checked by unit tests, so define them
ERR_CDEF_INCLASS = 'Cannot assign default value to fields in cdef classes, structs or unions'
ERR_BUF_DEFAULTS = 'Invalid buffer defaults specification (see docs)'
ERR_INVALID_SPECIALATTR_TYPE = 'Special attributes must not have a type declared'
class PostParse(ScopeTrackingTransform):
    """
    Basic interpretation of the parse tree, as well as validity
    checking that can be done on a very basic level on the parse
    tree (while still not being a problem with the basic syntax,
    as such).

    Specifically:
    - Default values to cdef assignments are turned into single
    assignments following the declaration (everywhere but in class
    bodies, where they raise a compile error)

    - Interpret some node structures into Python runtime values.
    Some nodes take compile-time arguments (currently:
    TemplatedTypeNode[args] and __cythonbufferdefaults__ = {args}),
    which should be interpreted. This happens in a general way
    and other steps should be taken to ensure validity.

    Type arguments cannot be interpreted in this way.

    - For __cythonbufferdefaults__ the arguments are checked for
    validity.

    TemplatedTypeNode has its directives interpreted:
    Any first positional argument goes into the "dtype" attribute,
    any "ndim" keyword argument goes into the "ndim" attribute and
    so on. Also it is checked that the directive combination is valid.
    - __cythonbufferdefaults__ attributes are parsed and put into the
    type information.

    Note: Currently Parsing.py does a lot of interpretation and
    reorganization that can be refactored into this transform
    if a more pure Abstract Syntax Tree is wanted.

    - Some invalid uses of := assignment expressions are detected
    """
    def __init__(self, context):
        super().__init__(context)
        self.specialattribute_handlers = {
            '__cythonbufferdefaults__' : self.handle_bufferdefaults
        }
        self.in_pattern_node = False

    def visit_LambdaNode(self, node):
        # unpack a lambda expression into the corresponding DefNode
        collector = YieldNodeCollector()
        collector.visitchildren(node.result_expr)
        if collector.has_yield or collector.has_await or isinstance(node.result_expr, ExprNodes.YieldExprNode):
            body = Nodes.ExprStatNode(
                node.result_expr.pos, expr=node.result_expr)
        else:
            body = Nodes.ReturnStatNode(
                node.result_expr.pos, value=node.result_expr)
        node.def_node = Nodes.DefNode(
            node.pos, name=node.name,
            args=node.args, star_arg=node.star_arg,
            starstar_arg=node.starstar_arg,
            body=body, doc=None)
        self.visitchildren(node)
        return node

    def visit_GeneratorExpressionNode(self, node):
        # unpack a generator expression into the corresponding DefNode
        collector = YieldNodeCollector()
        collector.visitchildren(node.loop, attrs=None, exclude=["iterator"])
        node.def_node = Nodes.DefNode(
            node.pos, name=node.name, doc=None,
            args=[], star_arg=None, starstar_arg=None,
            body=node.loop, is_async_def=collector.has_await,
            is_generator_expression=True)
        _AssignmentExpressionChecker.do_checks(node.loop, scope_is_class=self.scope_type in ("pyclass", "cclass"))
        self.visitchildren(node)
        return node

    def visit_ComprehensionNode(self, node):
        # enforce local scope also in Py2 for async generators (seriously, that's a Py3.6 feature...)
        if not node.has_local_scope:
            collector = YieldNodeCollector()
            collector.visitchildren(node.loop)
            if collector.has_await:
                node.has_local_scope = True
        _AssignmentExpressionChecker.do_checks(node.loop, scope_is_class=self.scope_type in ("pyclass", "cclass"))
        self.visitchildren(node)
        return node

    # cdef variables
    def handle_bufferdefaults(self, decl):
        if not isinstance(decl.default, ExprNodes.DictNode):
            raise PostParseError(decl.pos, ERR_BUF_DEFAULTS)
        self.scope_node.buffer_defaults_node = decl.default
        self.scope_node.buffer_defaults_pos = decl.pos

    def visit_CVarDefNode(self, node):
        # This assumes only plain names and pointers are assignable on
        # declaration. Also, it makes use of the fact that a cdef decl
        # must appear before the first use, so we don't have to deal with
        # "i = 3; cdef int i = i" and can simply move the nodes around.
        try:
            self.visitchildren(node)
            stats = [node]
            newdecls = []
            for decl in node.declarators:
                declbase = decl
                while isinstance(declbase, (Nodes.CPtrDeclaratorNode, Nodes.CConstDeclaratorNode)):
                    declbase = declbase.base
                if isinstance(declbase, Nodes.CNameDeclaratorNode):
                    if declbase.default is not None:
                        if self.scope_type in ('cclass', 'pyclass', 'struct'):
                            if isinstance(self.scope_node, Nodes.CClassDefNode):
                                handler = self.specialattribute_handlers.get(decl.name)
                                if handler:
                                    if decl is not declbase:
                                        raise PostParseError(decl.pos, ERR_INVALID_SPECIALATTR_TYPE)
                                    handler(decl)
                                    continue  # Remove declaration
                            raise PostParseError(decl.pos, ERR_CDEF_INCLASS)
                        first_assignment = self.scope_type != 'module'
                        stats.append(Nodes.SingleAssignmentNode(node.pos,
                            lhs=ExprNodes.NameNode(node.pos, name=declbase.name),
                            rhs=declbase.default, first=first_assignment,
                            from_pxd_cvardef=node.in_pxd))
                        declbase.default = None
                newdecls.append(decl)
            node.declarators = newdecls
            return stats
        except PostParseError as e:
            # An error in a cdef clause is ok, simply remove the declaration
            # and try to move on to report more errors
            self.context.nonfatal_error(e)
            return None

    # Split parallel assignments (a,b = b,a) into separate partial
    # assignments that are executed rhs-first using temps.  This
    # restructuring must be applied before type analysis so that known
    # types on rhs and lhs can be matched directly.  It is required in
    # the case that the types cannot be coerced to a Python type in
    # order to assign from a tuple.

    def visit_SingleAssignmentNode(self, node):
        self.visitchildren(node)
        return self._visit_assignment_node(node, [node.lhs, node.rhs])

    def visit_CascadedAssignmentNode(self, node):
        self.visitchildren(node)
        return self._visit_assignment_node(node, node.lhs_list + [node.rhs])

    def _visit_assignment_node(self, node, expr_list):
        """Flatten parallel assignments into separate single
        assignments or cascaded assignments.
        """
        if sum([ 1 for expr in expr_list
                 if expr.is_sequence_constructor or expr.is_string_literal ]) < 2:
            # no parallel assignments => nothing to do
            return node

        expr_list_list = []
        flatten_parallel_assignments(expr_list, expr_list_list)
        temp_refs = []
        eliminate_rhs_duplicates(expr_list_list, temp_refs)

        nodes = []
        for expr_list in expr_list_list:
            lhs_list = expr_list[:-1]
            rhs = expr_list[-1]
            if len(lhs_list) == 1:
                node = Nodes.SingleAssignmentNode(rhs.pos,
                    lhs = lhs_list[0], rhs = rhs)
            else:
                node = Nodes.CascadedAssignmentNode(rhs.pos,
                    lhs_list = lhs_list, rhs = rhs)
            nodes.append(node)

        if len(nodes) == 1:
            assign_node = nodes[0]
        else:
            assign_node = Nodes.ParallelAssignmentNode(nodes[0].pos, stats = nodes)

        if temp_refs:
            duplicates_and_temps = [ (temp.expression, temp)
                                     for temp in temp_refs ]
            sort_common_subsequences(duplicates_and_temps)
            for _, temp_ref in duplicates_and_temps[::-1]:
                assign_node = LetNode(temp_ref, assign_node)

        return assign_node

    def _flatten_sequence(self, seq, result):
        for arg in seq.args:
            if arg.is_sequence_constructor:
                self._flatten_sequence(arg, result)
            else:
                result.append(arg)
        return result

    def visit_DelStatNode(self, node):
        self.visitchildren(node)
        node.args = self._flatten_sequence(node, [])
        return node

    def visit_ExceptClauseNode(self, node):
        if node.is_except_as:
            # except-as must delete NameNode target at the end
            del_target = Nodes.DelStatNode(
                node.pos,
                args=[ExprNodes.NameNode(
                    node.target.pos, name=node.target.name)],
                ignore_nonexisting=True)
            node.body = Nodes.StatListNode(
                node.pos,
                stats=[Nodes.TryFinallyStatNode(
                    node.pos,
                    body=node.body,
                    finally_clause=Nodes.StatListNode(
                        node.pos,
                        stats=[del_target]))])
        self.visitchildren(node)
        return node

    def visit_AssertStatNode(self, node):
        """Extract the exception raising into a RaiseStatNode to simplify GIL handling.
        """
        if node.exception is None:
            node.exception = Nodes.RaiseStatNode(
                node.pos,
                exc_type=ExprNodes.NameNode(node.pos, name=EncodedString("AssertionError")),
                exc_value=node.value,
                exc_tb=None,
                cause=None,
                builtin_exc_name="AssertionError",
                wrap_tuple_value=True,
            )
            node.value = None
        self.visitchildren(node)
        return node

    def visit_ErrorNode(self, node):
        error(node.pos, node.what)
        return None

    def visit_MatchCaseNode(self, node):
        node.validate_targets()
        self.visitchildren(node)
        return node

    def visit_MatchNode(self, node):
        node.validate_irrefutable()
        self.visitchildren(node)
        return node

    def visit_PatternNode(self, node):
        in_pattern_node, self.in_pattern_node = self.in_pattern_node, True
        self.visitchildren(node)
        self.in_pattern_node = in_pattern_node
        return node

    def visit_JoinedStrNode(self, node):
        if self.in_pattern_node:
            error(node.pos, "f-strings are not accepted for pattern matching")
        self.visitchildren(node)
        return node

    def visit_DefNode(self, node):
        if (self.scope_type == "cclass" and
                node.name in ["__getreadbuffer__", "__getwritebuffer__", "__getsegcount__", "__getcharbuffer__"]):
            warning(node.pos, f"'{node.name}' relates to the old Python 2 buffer protocol "
                    "and is no longer used.", 2)
            return None  # drop the node - the arguments are invalid for a def node
        return self.visit_FuncDefNode(node)


class _AssignmentExpressionTargetNameFinder(TreeVisitor):
    def __init__(self):
        super().__init__()
        self.target_names = {}

    def find_target_names(self, target):
        if target.is_name:
            return [target.name]
        elif target.is_sequence_constructor:
            names = []
            for arg in target.args:
                names.extend(self.find_target_names(arg))
            return names
        # other targets are possible, but it isn't necessary to investigate them here
        return []

    def visit_ForInStatNode(self, node):
        self.target_names[node] = tuple(self.find_target_names(node.target))
        self.visitchildren(node)

    def visit_ComprehensionNode(self, node):
        pass  # don't recurse into nested comprehensions

    def visit_LambdaNode(self, node):
        pass  # don't recurse into nested lambdas/generator expressions

    def visit_Node(self, node):
        self.visitchildren(node)


class _AssignmentExpressionChecker(TreeVisitor):
    """
    Enforces rules on AssignmentExpressions within generator expressions and comprehensions
    """
    def __init__(self, loop_node, scope_is_class):
        super().__init__()

        target_name_finder = _AssignmentExpressionTargetNameFinder()
        target_name_finder.visit(loop_node)
        self.target_names_dict = target_name_finder.target_names
        self.in_iterator = False
        self.in_nested_generator = False
        self.scope_is_class = scope_is_class
        self.current_target_names = ()
        self.all_target_names = set()
        for names in self.target_names_dict.values():
            self.all_target_names.update(names)

    def _reset_state(self):
        old_state = (self.in_iterator, self.in_nested_generator, self.scope_is_class, self.all_target_names, self.current_target_names)
        # note: not resetting self.in_iterator here, see visit_LambdaNode() below
        self.in_nested_generator = False
        self.scope_is_class = False
        self.current_target_names = ()
        self.all_target_names = set()
        return old_state

    def _set_state(self, old_state):
        self.in_iterator, self.in_nested_generator, self.scope_is_class, self.all_target_names, self.current_target_names = old_state

    @classmethod
    def do_checks(cls, loop_node, scope_is_class):
        checker = cls(loop_node, scope_is_class)
        checker.visit(loop_node)

    def visit_ForInStatNode(self, node):
        if self.in_nested_generator:
            self.visitchildren(node)  # once nested, don't do anything special
            return

        current_target_names = self.current_target_names
        target_name = self.target_names_dict.get(node, None)
        if target_name:
            self.current_target_names += target_name

        self.in_iterator = True
        self.visit(node.iterator)
        self.in_iterator = False
        self.visitchildren(node, exclude=("iterator",))

        self.current_target_names = current_target_names

    def visit_AssignmentExpressionNode(self, node):
        if self.in_iterator:
            error(node.pos, "assignment expression cannot be used in a comprehension iterable expression")
        if self.scope_is_class:
            error(node.pos, "assignment expression within a comprehension cannot be used in a class body")
        if node.target_name in self.current_target_names:
            error(node.pos, "assignment expression cannot rebind comprehension iteration variable '%s'" %
                  node.target_name)
        elif node.target_name in self.all_target_names:
            error(node.pos, "comprehension inner loop cannot rebind assignment expression target '%s'" %
                  node.target_name)

    def visit_LambdaNode(self, node):
        # Don't reset "in_iterator" - an assignment expression in a lambda in an
        # iterator is explicitly tested by the Python testcases and banned.
        old_state = self._reset_state()
        # the lambda node's "def_node" is not set up at this point, so we need to recurse into it explicitly.
        self.visit(node.result_expr)
        self._set_state(old_state)

    def visit_ComprehensionNode(self, node):
        in_nested_generator = self.in_nested_generator
        self.in_nested_generator = True
        self.visitchildren(node)
        self.in_nested_generator = in_nested_generator

    def visit_GeneratorExpressionNode(self, node):
        in_nested_generator = self.in_nested_generator
        self.in_nested_generator = True
        # def_node isn't set up yet, so we need to visit the loop directly.
        self.visit(node.loop)
        self.in_nested_generator = in_nested_generator

    def visit_Node(self, node):
        self.visitchildren(node)


def eliminate_rhs_duplicates(expr_list_list, ref_node_sequence):
    """Replace rhs items by LetRefNodes if they appear more than once.
    Creates a sequence of LetRefNodes that set up the required temps
    and appends them to ref_node_sequence.  The input list is modified
    in-place.
    """
    seen_nodes = set()
    ref_nodes = {}
    def find_duplicates(node):
        if node.is_literal or node.is_name:
            # no need to replace those; can't include attributes here
            # as their access is not necessarily side-effect free
            return
        if node in seen_nodes:
            if node not in ref_nodes:
                ref_node = LetRefNode(node)
                ref_nodes[node] = ref_node
                ref_node_sequence.append(ref_node)
        else:
            seen_nodes.add(node)
            if node.is_sequence_constructor:
                for item in node.args:
                    find_duplicates(item)

    for expr_list in expr_list_list:
        rhs = expr_list[-1]
        find_duplicates(rhs)
    if not ref_nodes:
        return

    def substitute_nodes(node):
        if node in ref_nodes:
            return ref_nodes[node]
        elif node.is_sequence_constructor:
            node.args = list(map(substitute_nodes, node.args))
        return node

    # replace nodes inside of the common subexpressions
    for node in ref_nodes:
        if node.is_sequence_constructor:
            node.args = list(map(substitute_nodes, node.args))

    # replace common subexpressions on all rhs items
    for expr_list in expr_list_list:
        expr_list[-1] = substitute_nodes(expr_list[-1])

def sort_common_subsequences(items):
    """Sort items/subsequences so that all items and subsequences that
    an item contains appear before the item itself.  This is needed
    because each rhs item must only be evaluated once, so its value
    must be evaluated first and then reused when packing sequences
    that contain it.

    This implies a partial order, and the sort must be stable to
    preserve the original order as much as possible, so we use a
    simple insertion sort (which is very fast for short sequences, the
    normal case in practice).
    """
    def contains(seq, x):
        for item in seq:
            if item is x:
                return True
            elif item.is_sequence_constructor and contains(item.args, x):
                return True
        return False
    def lower_than(a,b):
        return b.is_sequence_constructor and contains(b.args, a)

    for pos, item in enumerate(items):
        key = item[1]  # the ResultRefNode which has already been injected into the sequences
        new_pos = pos
        for i in range(pos-1, -1, -1):
            if lower_than(key, items[i][0]):
                new_pos = i
        if new_pos != pos:
            for i in range(pos, new_pos, -1):
                items[i] = items[i-1]
            items[new_pos] = item


def unpack_string_to_character_literals(literal):
    chars = []
    pos = literal.pos
    stype = literal.__class__
    sval = literal.value
    sval_type = sval.__class__
    for char in sval:
        cval = sval_type(char)
        chars.append(stype(pos, value=cval))
    return chars


@cython.cfunc
def flatten_parallel_assignments(input: list, output: list):
    #  The input is a list of expression nodes, representing the LHSs
    #  and RHS of one (possibly cascaded) assignment statement.  For
    #  sequence constructors, rearranges the matching parts of both
    #  sides into a list of equivalent assignments between the
    #  individual elements.  This transformation is applied
    #  recursively, so that nested structures get matched as well.
    rhs = input[-1]
    if (not (rhs.is_sequence_constructor or isinstance(rhs, ExprNodes.UnicodeNode))
            or not sum([lhs.is_sequence_constructor for lhs in input[:-1]])):
        output.append(input)
        return

    complete_assignments = []

    if rhs.is_sequence_constructor:
        rhs_args = rhs.args
    elif rhs.is_string_literal:
        rhs_args = unpack_string_to_character_literals(rhs)

    starred_targets: cython.Py_ssize_t
    lhs_size: cython.Py_ssize_t
    rhs_size: cython.Py_ssize_t = len(rhs_args)
    lhs_targets = [[] for _ in range(rhs_size)]
    starred_assignments = []

    for lhs in input[:-1]:
        if not lhs.is_sequence_constructor:
            if lhs.is_starred:
                error(lhs.pos, "starred assignment target must be in a list or tuple")
            complete_assignments.append(lhs)
            continue
        lhs_size = len(lhs.args)
        starred_targets = 0
        for expr in lhs.args:
            starred_targets += bool(expr.is_starred)
        if starred_targets > 1:
            error(lhs.pos, "more than 1 starred expression in assignment")
            output.append([lhs,rhs])
            continue
        elif lhs_size - starred_targets > rhs_size:
            error(lhs.pos, "need more than %d value%s to unpack"
                  % (rhs_size, (rhs_size != 1) and 's' or ''))
            output.append([lhs,rhs])
            continue
        elif starred_targets:
            map_starred_assignment(lhs_targets, starred_assignments,
                                   lhs.args, rhs_args)
        elif lhs_size < rhs_size:
            error(lhs.pos, "too many values to unpack (expected %d, got %d)"
                  % (lhs_size, rhs_size))
            output.append([lhs,rhs])
            continue
        else:
            for targets, expr in zip(lhs_targets, lhs.args):
                targets.append(expr)

    if complete_assignments:
        complete_assignments.append(rhs)
        output.append(complete_assignments)

    # recursively flatten partial assignments
    for cascade, rhs in zip(lhs_targets, rhs_args):
        if cascade:
            cascade.append(rhs)
            flatten_parallel_assignments(cascade, output)

    # recursively flatten starred assignments
    for cascade in starred_assignments:
        if cascade[0].is_sequence_constructor:
            flatten_parallel_assignments(cascade, output)
        else:
            output.append(cascade)


@cython.cfunc
def map_starred_assignment(lhs_targets: list, starred_assignments: list, lhs_args: list, rhs_args: list):
    # Appends the fixed-position LHS targets to the target list that
    # appear left and right of the starred argument.
    #
    # The starred_assignments list receives a new tuple
    # (lhs_target, rhs_values_list) that maps the remaining arguments
    # (those that match the starred target) to a list.

    # left side of the starred target
    i: cython.Py_ssize_t
    starred: cython.Py_ssize_t
    lhs_remaining: cython.Py_ssize_t
    for i, (targets, expr) in enumerate(zip(lhs_targets, lhs_args)):
        if expr.is_starred:
            starred = i
            lhs_remaining = len(lhs_args) - i - 1
            break
        targets.append(expr)
    else:
        raise InternalError("no starred arg found when splitting starred assignment")

    # right side of the starred target
    for i, (targets, expr) in enumerate(zip(lhs_targets[-lhs_remaining:],
                                            lhs_args[starred + 1:])):
        targets.append(expr)

    # the starred target itself, must be assigned a (potentially empty) list
    target = lhs_args[starred].target  # unpack starred node
    starred_rhs = rhs_args[starred:]
    if lhs_remaining:
        starred_rhs = starred_rhs[:-lhs_remaining]
    if starred_rhs:
        pos = starred_rhs[0].pos
    else:
        pos = target.pos
    starred_assignments.append([
        target, ExprNodes.ListNode(pos=pos, args=starred_rhs)])


class PxdPostParse(CythonTransform, SkipDeclarations):
    """
    Basic interpretation/validity checking that should only be
    done on pxd trees.

    A lot of this checking currently happens in the parser; but
    what is listed below happens here.

    - "def" functions are let through only if they fill the
    getbuffer/releasebuffer slots

    - cdef functions are let through only if they are on the
    top level and are declared "inline"
    """
    ERR_INLINE_ONLY = "function definition in pxd file must be declared 'cdef inline'"
    ERR_NOGO_WITH_INLINE = "inline function definition in pxd file cannot be '%s'"

    def __call__(self, node):
        self.scope_type = 'pxd'
        return super().__call__(node)

    def visit_CClassDefNode(self, node):
        old = self.scope_type
        self.scope_type = 'cclass'
        self.visitchildren(node)
        self.scope_type = old
        return node

    def visit_FuncDefNode(self, node):
        # FuncDefNode always come with an implementation (without
        # an imp they are CVarDefNodes..)
        err = self.ERR_INLINE_ONLY

        if (isinstance(node, Nodes.DefNode) and self.scope_type == 'cclass'
                and node.name in ('__getbuffer__', '__releasebuffer__')):
            err = None  # allow these slots

        if isinstance(node, Nodes.CFuncDefNode):
            if ('inline' in node.modifiers and
                    self.scope_type in ('pxd', 'cclass')):
                node.inline_in_pxd = True
                if node.visibility != 'private':
                    err = self.ERR_NOGO_WITH_INLINE % node.visibility
                elif node.api:
                    err = self.ERR_NOGO_WITH_INLINE % 'api'
                else:
                    err = None  # allow inline function
            else:
                err = self.ERR_INLINE_ONLY

        if err:
            self.context.nonfatal_error(PostParseError(node.pos, err))
            return None
        else:
            return node


class TrackNumpyAttributes(VisitorTransform, SkipDeclarations):
    # TODO: Make name handling as good as in InterpretCompilerDirectives() below - probably best to merge the two.
    def __init__(self):
        super().__init__()
        self.numpy_module_names = set()

    def visit_CImportStatNode(self, node):
        if node.module_name == "numpy":
            self.numpy_module_names.add(node.as_name or "numpy")
        return node

    def visit_AttributeNode(self, node):
        self.visitchildren(node)
        obj = node.obj
        if (obj.is_name and obj.name in self.numpy_module_names) or obj.is_numpy_attribute:
            node.is_numpy_attribute = True
        return node

    visit_Node = VisitorTransform.recurse_to_children


class InterpretCompilerDirectives(CythonTransform):
    """
    After parsing, directives can be stored in a number of places:
    - #cython-comments at the top of the file (stored in ModuleNode)
    - Command-line arguments overriding these
    - @cython.directivename decorators
    - with cython.directivename: statements
    - replaces "cython.compiled" with BoolNode(value=True)
      allowing unreachable blocks to be removed at a fairly early stage
      before cython typing rules are forced on applied

    This transform is responsible for interpreting these various sources
    and store the directive in two ways:
    - Set the directives attribute of the ModuleNode for global directives.
    - Use a CompilerDirectivesNode to override directives for a subtree.

    (The first one is primarily to not have to modify with the tree
    structure, so that ModuleNode stay on top.)

    The directives are stored in dictionaries from name to value in effect.
    Each such dictionary is always filled in for all possible directives,
    using default values where no value is given by the user.

    The available directives are controlled in Options.py.

    Note that we have to run this prior to analysis, and so some minor
    duplication of functionality has to occur: We manually track cimports
    and which names the "cython" module may have been imported to.
    """
    unop_method_nodes = {
        'typeof': ExprNodes.TypeofNode,

        'operator.address': ExprNodes.AmpersandNode,
        'operator.dereference': ExprNodes.DereferenceNode,
        'operator.preincrement' : ExprNodes.inc_dec_constructor(True, '++'),
        'operator.predecrement' : ExprNodes.inc_dec_constructor(True, '--'),
        'operator.postincrement': ExprNodes.inc_dec_constructor(False, '++'),
        'operator.postdecrement': ExprNodes.inc_dec_constructor(False, '--'),
        'operator.typeid'       : ExprNodes.TypeidNode,

        # For backwards compatibility.
        'address': ExprNodes.AmpersandNode,
    }

    binop_method_nodes = {
        'operator.comma'        : ExprNodes.c_binop_constructor(','),
    }

    special_methods = {
        'declare', 'union', 'struct', 'typedef',
        'sizeof', 'cast', 'pointer', 'compiled',
        'NULL', 'fused_type', 'parallel',
    }
    special_methods.update(unop_method_nodes)

    valid_cython_submodules = {
        'cimports',
        'dataclasses',
        'operator',
        'parallel',
        'view',
    }

    valid_parallel_directives = {
        "parallel",
        "prange",
        "threadid",
        #"threadsavailable",
    }

    def __init__(self, context, compilation_directive_defaults):
        super().__init__(context)
        self.cython_module_names = set()
        self.directive_names = {'staticmethod': 'staticmethod'}
        self.parallel_directives = {}
        directives = copy.deepcopy(Options.get_directive_defaults())
        for key, value in compilation_directive_defaults.items():
            directives[str(key)] = copy.deepcopy(value)
        self.directives = directives

    def check_directive_scope(self, pos, directive, scope):
        legal_scopes = Options.directive_scopes.get(directive, None)
        if legal_scopes and scope not in legal_scopes:
            self.context.nonfatal_error(PostParseError(pos, 'The %s compiler directive '
                                        'is not allowed in %s scope' % (directive, scope)))
            return False
        else:
            if directive not in Options.directive_types:
                error(pos, "Invalid directive: '%s'." % (directive,))
            return True

    def _check_valid_cython_module(self, pos, module_name):
        if not module_name.startswith("cython."):
            return
        submodule = module_name.split('.', 2)[1]
        if submodule in self.valid_cython_submodules:
            return

        extra = ""
        # This is very rarely used, so don't waste space on static tuples.
        hints = [
            line.split() for line in """\
                imp                  cimports
                cimp                 cimports
                para                 parallel
                parra                parallel
                dataclass            dataclasses
            """.splitlines()[:-1]
        ]
        for wrong, correct in hints:
            if module_name.startswith("cython." + wrong):
                extra = "Did you mean 'cython.%s' ?" % correct
                break
        if not extra:
            is_simple_cython_name = submodule in Options.directive_types
            if not is_simple_cython_name and not submodule.startswith("_"):
                # Try to find it in the Shadow module (i.e. the pure Python namespace of cython.*).
                # FIXME: use an internal reference of "cython.*" names instead of Shadow.py
                from .. import Shadow
                is_simple_cython_name = hasattr(Shadow, submodule)
            if is_simple_cython_name:
                extra = "Instead, use 'import cython' and then 'cython.%s'." % submodule

        error(pos, "'%s' is not a valid cython.* module%s%s" % (
            module_name,
            ". " if extra else "",
            extra,
        ))

    # Set up processing and handle the cython: comments.
    def visit_ModuleNode(self, node):
        for key in sorted(node.directive_comments):
            if not self.check_directive_scope(node.pos, key, 'module'):
                self.wrong_scope_error(node.pos, key, 'module')
                del node.directive_comments[key]

        self.module_scope = node.scope

        self.directives.update(node.directive_comments)
        node.directives = self.directives
        node.parallel_directives = self.parallel_directives
        self.visitchildren(node)
        node.cython_module_names = self.cython_module_names
        return node

    def visit_CompilerDirectivesNode(self, node):
        old_directives, self.directives = self.directives, node.directives
        self.visitchildren(node)
        self.directives = old_directives
        return node

    # The following four functions track imports and cimports that
    # begin with "cython"
    def is_cython_directive(self, name):
        return (name in Options.directive_types or
                name in self.special_methods or
                PyrexTypes.parse_basic_type(name))

    def is_parallel_directive(self, full_name, pos):
        """
        Checks to see if fullname (e.g. cython.parallel.prange) is a valid
        parallel directive. If it is a star import it also updates the
        parallel_directives.
        """
        result = (full_name + ".").startswith("cython.parallel.")

        if result:
            directive = full_name.split('.')
            if full_name == "cython.parallel":
                self.parallel_directives["parallel"] = "cython.parallel"
            elif full_name == "cython.parallel.*":
                for name in self.valid_parallel_directives:
                    self.parallel_directives[name] = "cython.parallel.%s" % name
            elif (len(directive) != 3 or
                  directive[-1] not in self.valid_parallel_directives):
                error(pos, "No such directive: %s" % full_name)

        return result

    def visit_CImportStatNode(self, node):
        module_name = node.module_name
        if module_name == "cython.cimports":
            error(node.pos, "Cannot cimport the 'cython.cimports' package directly, only submodules.")
        if module_name.startswith("cython.cimports."):
            if node.as_name and node.as_name != 'cython':
                node.module_name = module_name[len("cython.cimports."):]
                return node
            error(node.pos,
                  "Python cimports must use 'from cython.cimports... import ...'"
                  " or 'import ... as ...', not just 'import ...'")

        if module_name == "cython":
            self.cython_module_names.add(node.as_name or "cython")
        elif module_name.startswith("cython."):
            if module_name.startswith("cython.parallel."):
                error(node.pos, node.module_name + " is not a module")
            else:
                self._check_valid_cython_module(node.pos, module_name)

            if module_name == "cython.parallel":
                if node.as_name and node.as_name != "cython":
                    self.parallel_directives[node.as_name] = module_name
                else:
                    self.cython_module_names.add("cython")
                    self.parallel_directives[
                                    "cython.parallel"] = module_name
            elif node.as_name:
                self.directive_names[node.as_name] = module_name[7:]
            else:
                self.cython_module_names.add("cython")
            # if this cimport was a compiler directive, we don't
            # want to leave the cimport node sitting in the tree
            return None
        return node

    def visit_FromCImportStatNode(self, node):
        module_name = node.module_name
        if module_name == "cython.cimports" or module_name.startswith("cython.cimports."):
            # only supported for convenience
            return self._create_cimport_from_import(
                node.pos, module_name, node.relative_level, node.imported_names)
        elif not node.relative_level and (
                module_name == "cython" or module_name.startswith("cython.")):
            self._check_valid_cython_module(node.pos, module_name)
            submodule = (module_name + ".")[7:]
            newimp = []
            for pos, name, as_name in node.imported_names:
                full_name = submodule + name
                qualified_name = "cython." + full_name
                if self.is_parallel_directive(qualified_name, node.pos):
                    # from cython cimport parallel, or
                    # from cython.parallel cimport parallel, prange, ...
                    self.parallel_directives[as_name or name] = qualified_name
                elif self.is_cython_directive(full_name):
                    self.directive_names[as_name or name] = full_name
                elif full_name in ['dataclasses', 'typing']:
                    self.directive_names[as_name or name] = full_name
                    # unlike many directives, still treat it as a regular module
                    newimp.append((pos, name, as_name))
                else:
                    newimp.append((pos, name, as_name))

            if not newimp:
                return None

            node.imported_names = newimp
        return node

    def visit_FromImportStatNode(self, node):
        import_node = node.module
        module_name = import_node.module_name.value
        if module_name == "cython.cimports" or module_name.startswith("cython.cimports."):
            imported_names = []
            for name, name_node in node.items:
                imported_names.append(
                    (name_node.pos, name, None if name == name_node.name else name_node.name))
            return self._create_cimport_from_import(
                node.pos, module_name, import_node.level, imported_names)
        elif module_name == "cython" or module_name.startswith("cython."):
            self._check_valid_cython_module(import_node.module_name.pos, module_name)
            submodule = (module_name + ".")[7:]
            newimp = []
            for name, name_node in node.items:
                full_name = submodule + name
                qualified_name = "cython." + full_name
                if self.is_parallel_directive(qualified_name, node.pos):
                    self.parallel_directives[name_node.name] = qualified_name
                elif self.is_cython_directive(full_name):
                    self.directive_names[name_node.name] = full_name
                else:
                    newimp.append((name, name_node))
            if not newimp:
                return None
            node.items = newimp
        return node

    def _create_cimport_from_import(self, node_pos, module_name, level, imported_names):
        if module_name == "cython.cimports" or module_name.startswith("cython.cimports."):
            module_name = EncodedString(module_name[len("cython.cimports."):])  # may be empty

        if module_name:
            # from cython.cimports.a.b import x, y, z  =>  from a.b cimport x, y, z
            return Nodes.FromCImportStatNode(
                node_pos, module_name=module_name,
                relative_level=level,
                imported_names=imported_names)
        else:
            # from cython.cimports import x, y, z  =>  cimport x; cimport y; cimport z
            return [
                Nodes.CImportStatNode(
                    pos,
                    module_name=dotted_name,
                    as_name=as_name,
                    is_absolute=level == 0)
                for pos, dotted_name, as_name in imported_names
            ]

    def visit_SingleAssignmentNode(self, node):
        if isinstance(node.rhs, ExprNodes.ImportNode):
            module_name = node.rhs.module_name.value
            if module_name != "cython" and not module_name.startswith("cython."):
                return node

            node = Nodes.CImportStatNode(node.pos, module_name=module_name, as_name=node.lhs.name)
            node = self.visit_CImportStatNode(node)
        else:
            self.visitchildren(node)

        return node

    def visit_NameNode(self, node):
        if node.annotation:
            self.visitchild(node, 'annotation')
        if node.name in self.cython_module_names:
            node.is_cython_module = True
        else:
            directive = self.directive_names.get(node.name)
            if directive is not None:
                node.cython_attribute = directive
        if node.as_cython_attribute() == "compiled":
            return ExprNodes.BoolNode(node.pos, value=True)  # replace early so unused branches can be dropped
                # before they have a chance to cause compile-errors
        return node

    def visit_AttributeNode(self, node):
        self.visitchildren(node)
        if node.as_cython_attribute() == "compiled":
            return ExprNodes.BoolNode(node.pos, value=True)  # replace early so unused branches can be dropped
                # before they have a chance to cause compile-errors
        return node

    def visit_AnnotationNode(self, node):
        # for most transforms annotations are left unvisited (because they're unevaluated)
        # however, it is important to pick up compiler directives from them
        if node.expr:
            self.visit(node.expr)
        return node

    def visit_NewExprNode(self, node):
        self.visitchild(node, 'cppclass')
        self.visitchildren(node)
        return node

    def try_to_parse_directives(self, node):
        # If node is the contents of an directive (in a with statement or
        # decorator), returns a list of (directivename, value) pairs.
        # Otherwise, returns None
        if isinstance(node, ExprNodes.CallNode):
            self.visitchild(node, 'function')
            optname = node.function.as_cython_attribute()
            if optname:
                directivetype = Options.directive_types.get(optname)
                if directivetype:
                    args, kwds = node.explicit_args_kwds()
                    directives = []
                    key_value_pairs = []
                    if kwds is not None and directivetype is not dict:
                        for keyvalue in kwds.key_value_pairs:
                            key, value = keyvalue
                            sub_optname = "%s.%s" % (optname, key.value)
                            if Options.directive_types.get(sub_optname):
                                directives.append(self.try_to_parse_directive(sub_optname, [value], None, keyvalue.pos))
                            else:
                                key_value_pairs.append(keyvalue)
                        if not key_value_pairs:
                            kwds = None
                        else:
                            kwds.key_value_pairs = key_value_pairs
                        if directives and not kwds and not args:
                            return directives
                    directives.append(self.try_to_parse_directive(optname, args, kwds, node.function.pos))
                    return directives
        elif isinstance(node, (ExprNodes.AttributeNode, ExprNodes.NameNode)):
            self.visit(node)
            optname = node.as_cython_attribute()
            if optname:
                directivetype = Options.directive_types.get(optname)
                if directivetype is bool:
                    arg = ExprNodes.BoolNode(node.pos, value=True)
                    return [self.try_to_parse_directive(optname, [arg], None, node.pos)]
                elif directivetype is None or directivetype is Options.DEFER_ANALYSIS_OF_ARGUMENTS:
                    return [(optname, None)]
                else:
                    raise PostParseError(
                        node.pos, "The '%s' directive should be used as a function call." % optname)
        return None

    def try_to_parse_directive(self, optname, args, kwds, pos):
        if optname == 'np_pythran' and not self.context.cpp:
            raise PostParseError(pos, 'The %s directive can only be used in C++ mode.' % optname)
        elif optname == 'exceptval':
            # default: exceptval(None, check=True)
            arg_error = len(args) > 1
            check = True
            if kwds and kwds.key_value_pairs:
                kw = kwds.key_value_pairs[0]
                if (len(kwds.key_value_pairs) == 1 and
                        kw.key.is_string_literal and kw.key.value == 'check' and
                        isinstance(kw.value, ExprNodes.BoolNode)):
                    check = kw.value.value
                else:
                    arg_error = True
            if arg_error:
                raise PostParseError(
                    pos, 'The exceptval directive takes 0 or 1 positional arguments and the boolean keyword "check"')
            return ('exceptval', (args[0] if args else None, check))

        directivetype = Options.directive_types.get(optname)
        if len(args) == 1 and isinstance(args[0], ExprNodes.NoneNode):
            return optname, Options.get_directive_defaults()[optname]
        elif directivetype is bool:
            if kwds is not None or len(args) != 1 or not isinstance(args[0], ExprNodes.BoolNode):
                raise PostParseError(pos,
                    'The %s directive takes one compile-time boolean argument' % optname)
            return (optname, args[0].value)
        elif directivetype is int:
            if kwds is not None or len(args) != 1 or not isinstance(args[0], ExprNodes.IntNode):
                raise PostParseError(pos,
                    'The %s directive takes one compile-time integer argument' % optname)
            return (optname, int(args[0].value))
        elif directivetype is str:
            if kwds is not None or len(args) != 1 or not isinstance(args[0], ExprNodes.UnicodeNode):
                raise PostParseError(pos,
                    'The %s directive takes one compile-time string argument' % optname)
            return (optname, str(args[0].value))
        elif directivetype is type:
            if kwds is not None or len(args) != 1:
                raise PostParseError(pos,
                    'The %s directive takes one type argument' % optname)
            return (optname, args[0])
        elif directivetype is dict:
            if len(args) != 0:
                raise PostParseError(pos,
                    'The %s directive takes no prepositional arguments' % optname)
            return optname, kwds.as_python_dict()
        elif directivetype is list:
            if kwds and len(kwds.key_value_pairs) != 0:
                raise PostParseError(pos,
                    'The %s directive takes no keyword arguments' % optname)
            return optname, [ str(arg.value) for arg in args ]
        elif callable(directivetype):
            if kwds is not None or len(args) != 1 or not isinstance(args[0], ExprNodes.UnicodeNode):
                raise PostParseError(pos,
                    'The %s directive takes one compile-time string argument' % optname)
            return (optname, directivetype(optname, str(args[0].value)))
        elif directivetype is Options.DEFER_ANALYSIS_OF_ARGUMENTS:
            # signal to pass things on without processing
            return (optname, (args, kwds.as_python_dict() if kwds else {}))
        else:
            assert False

    def visit_with_directives(self, node, directives, contents_directives):
        # contents_directives may be None
        if not directives:
            assert not contents_directives
            return self.visit_Node(node)

        old_directives = self.directives
        new_directives = Options.copy_inherited_directives(old_directives, **directives)
        if contents_directives is not None:
            new_contents_directives = Options.copy_inherited_directives(
                old_directives, **contents_directives)
        else:
            new_contents_directives = new_directives

        if new_directives == old_directives:
            return self.visit_Node(node)

        self.directives = new_directives
        if (contents_directives is not None and
                new_contents_directives != new_directives):
            # we need to wrap the node body in a compiler directives node
            node.body = Nodes.StatListNode(
                node.body.pos,
                stats=[
                    Nodes.CompilerDirectivesNode(
                        node.body.pos,
                        directives=new_contents_directives,
                        body=node.body)
                ]
            )
        retbody = self.visit_Node(node)
        self.directives = old_directives

        if isinstance(retbody, Nodes.CompilerDirectivesNode):
            new_directives.update(retbody.directives)
            retbody = retbody.body
        if not isinstance(retbody, Nodes.StatListNode):
            retbody = Nodes.StatListNode(node.pos, stats=[retbody])
        return Nodes.CompilerDirectivesNode(
            retbody.pos, body=retbody, directives=new_directives, is_terminator=retbody.is_terminator)

    # Handle decorators
    def visit_FuncDefNode(self, node):
        directives, contents_directives = self._extract_directives(node, 'function')
        return self.visit_with_directives(node, directives, contents_directives)

    def visit_CVarDefNode(self, node):
        directives, _ = self._extract_directives(node, 'function')
        for name, value in directives.items():
            if name == 'locals':
                node.directive_locals = value
            elif name not in ('final', 'staticmethod'):
                self.context.nonfatal_error(PostParseError(
                    node.pos,
                    "Cdef functions can only take cython.locals(), "
                    "staticmethod, or final decorators, got %s." % name))
        return self.visit_with_directives(node, directives, contents_directives=None)

    def visit_CClassDefNode(self, node):
        directives, contents_directives = self._extract_directives(node, 'cclass')
        return self.visit_with_directives(node, directives, contents_directives)

    def visit_CppClassNode(self, node):
        directives, contents_directives = self._extract_directives(node, 'cppclass')
        return self.visit_with_directives(node, directives, contents_directives)

    def visit_PyClassDefNode(self, node):
        directives, contents_directives = self._extract_directives(node, 'class')
        return self.visit_with_directives(node, directives, contents_directives)

    def _extract_directives(self, node, scope_name):
        """
        Returns two dicts - directives applied to this function/class
        and directives applied to its contents. They aren't always the
        same (since e.g. cfunc should not be applied to inner functions)
        """
        if not node.decorators:
            return {}, {}
        # Split the decorators into two lists -- real decorators and directives
        directives = []
        realdecs = []
        both = []
        current_opt_dict = dict(self.directives)
        missing = object()
        # Decorators coming first take precedence.
        for dec in node.decorators[::-1]:
            new_directives = self.try_to_parse_directives(dec.decorator)
            if new_directives is not None:
                for directive in new_directives:
                    if self.check_directive_scope(node.pos, directive[0], scope_name):
                        name, value = directive
                        if name in ('nogil', 'with_gil'):
                            if value is None:
                                value = True
                            else:
                                args, kwds = value
                                if kwds or len(args) != 1 or not isinstance(args[0], ExprNodes.BoolNode):
                                    raise PostParseError(dec.pos, 'The %s directive takes one compile-time boolean argument' % name)
                                value = args[0].value
                            directive = (name, value)
                        if current_opt_dict.get(name, missing) != value:
                            if name == 'cfunc' and 'ufunc' in current_opt_dict:
                                error(dec.pos, "Cannot apply @cfunc to @ufunc, please reverse the decorators.")
                            directives.append(directive)
                            current_opt_dict[name] = value
                        else:
                            warning(dec.pos, "Directive does not change previous value (%s%s)" % (
                                name, '=%r' % value if value is not None else ''))
                        if directive[0] == 'staticmethod':
                            both.append(dec)
                    # Adapt scope type based on decorators that change it.
                    if directive[0] == 'cclass' and scope_name == 'class':
                        scope_name = 'cclass'
            else:
                realdecs.append(dec)
        node.decorators = realdecs[::-1] + both[::-1]
        # merge or override repeated directives
        optdict = {}
        contents_optdict = {}
        for name, value in directives:
            if name in optdict:
                old_value = optdict[name]
                # keywords and arg lists can be merged, everything
                # else overrides completely
                if isinstance(old_value, dict):
                    old_value.update(value)
                elif isinstance(old_value, list):
                    old_value.extend(value)
                else:
                    if name == "collection_type" and value != optdict[name]:
                        error(node.pos, "Multiple values of collection_type are not supported")
                    optdict[name] = value
            else:
                optdict[name] = value
            if name not in Options.immediate_decorator_directives:
                contents_optdict[name] = value
        return optdict, contents_optdict

    # Handle with-statements
    def visit_WithStatNode(self, node):
        directive_dict = {}
        for directive in self.try_to_parse_directives(node.manager) or []:
            if directive is None:
                continue
            if node.target is not None:
                self.context.nonfatal_error(
                    PostParseError(node.pos, "Compiler directive with statements cannot contain 'as'"))
                continue
            name, value = directive
            if name in ('nogil', 'gil'):
                # special case: in pure mode, "with nogil" spells "with cython.nogil"
                return self._transform_with_gil(node, name)
            elif name == "critical_section":
                args, kwds = value
                return self._transform_critical_section(node, args, kwds)
            elif self.check_directive_scope(node.pos, name, 'with statement'):
                directive_dict[name] = value
        if directive_dict:
            return self.visit_with_directives(node.body, directive_dict, contents_directives=None)
        return self.visit_Node(node)

    def _transform_with_gil(self, node, state):
        assert state in ('gil', 'nogil')
        manager = node.manager
        condition = None
        if isinstance(manager, ExprNodes.SimpleCallNode) and manager.args:
            if len(manager.args) > 1:
                self.context.nonfatal_error(
                    PostParseError(node.pos, "Compiler directive %s accepts one positional argument." % state))
            condition = manager.args[0]
        elif isinstance(manager, ExprNodes.GeneralCallNode):
            self.context.nonfatal_error(
                PostParseError(node.pos, "Compiler directive %s accepts one positional argument." % state))
        node = Nodes.GILStatNode(node.pos, state=state, body=node.body, condition=condition)
        return self.visit_Node(node)

    def _transform_critical_section(self, node, args, kwds):
        if len(args) < 1 or len(args) > 2 or kwds:
            self.context.nonfatal_error(
                PostParseError(node.pos, "critical_section directive accepts one or two positional arguments")
            )
        node = Nodes.CriticalSectionStatNode(
            node.pos, args=args, body=node.body
        )
        return self.visit_Node(node)


class ParallelRangeTransform(CythonTransform, SkipDeclarations):
    """
    Transform cython.parallel stuff. The parallel_directives come from the
    module node, set there by InterpretCompilerDirectives.

        x = cython.parallel.threadavailable()   -> ParallelThreadAvailableNode
        with nogil, cython.parallel.parallel(): -> ParallelWithBlockNode
            print cython.parallel.threadid()    -> ParallelThreadIdNode
            for i in cython.parallel.prange(...):  -> ParallelRangeNode
                ...
    """

    # a list of names, maps 'cython.parallel.prange' in the code to
    # ['cython', 'parallel', 'prange']
    parallel_directive = None

    # Indicates whether a namenode in an expression is the cython module
    namenode_is_cython_module = False

    # Keep track of whether we are the context manager of a 'with' statement
    in_context_manager_section = False

    # One of 'prange' or 'with parallel'. This is used to disallow closely
    # nested 'with parallel:' blocks
    state = None

    directive_to_node = {
        "cython.parallel.parallel": Nodes.ParallelWithBlockNode,
        # u"cython.parallel.threadsavailable": ExprNodes.ParallelThreadsAvailableNode,
        "cython.parallel.threadid": ExprNodes.ParallelThreadIdNode,
        "cython.parallel.prange": Nodes.ParallelRangeNode,
    }

    def node_is_parallel_directive(self, node):
        return node.name in self.parallel_directives or node.is_cython_module

    def get_directive_class_node(self, node):
        """
        Figure out which parallel directive was used and return the associated
        Node class.

        E.g. for a cython.parallel.prange() call we return ParallelRangeNode
        """
        if self.namenode_is_cython_module:
            directive = '.'.join(self.parallel_directive)
        else:
            directive = self.parallel_directives[self.parallel_directive[0]]
            directive = '%s.%s' % (directive,
                                   '.'.join(self.parallel_directive[1:]))
            directive = directive.rstrip('.')

        cls = self.directive_to_node.get(directive)
        if cls is None and not (self.namenode_is_cython_module and
                                self.parallel_directive[0] != 'parallel'):
            error(node.pos, "Invalid directive: %s" % directive)

        self.namenode_is_cython_module = False
        self.parallel_directive = None

        return cls

    def visit_ModuleNode(self, node):
        """
        If any parallel directives were imported, copy them over and visit
        the AST
        """
        if node.parallel_directives:
            self.parallel_directives = node.parallel_directives
            return self.visit_Node(node)

        # No parallel directives were imported, so they can't be used :)
        return node

    def visit_NameNode(self, node):
        if self.node_is_parallel_directive(node):
            self.parallel_directive = [node.name]
            self.namenode_is_cython_module = node.is_cython_module
        return node

    def visit_AttributeNode(self, node):
        self.visitchildren(node)
        if self.parallel_directive:
            self.parallel_directive.append(node.attribute)
        return node

    def visit_CallNode(self, node):
        self.visitchild(node, 'function')
        if not self.parallel_directive:
            self.visitchildren(node, exclude=('function',))
            return node

        # We are a parallel directive, replace this node with the
        # corresponding ParallelSomethingSomething node

        if isinstance(node, ExprNodes.GeneralCallNode):
            args = node.positional_args.args
            kwargs = node.keyword_args
        else:
            args = node.args
            kwargs = {}

        parallel_directive_class = self.get_directive_class_node(node)
        if parallel_directive_class:
            # Note: in case of a parallel() the body is set by
            # visit_WithStatNode
            node = parallel_directive_class(node.pos, args=args, kwargs=kwargs)

        return node

    def visit_WithStatNode(self, node):
        "Rewrite with cython.parallel.parallel() blocks"
        newnode = self.visit(node.manager)

        if isinstance(newnode, Nodes.ParallelWithBlockNode):
            if self.state == 'parallel with':
                error(node.manager.pos,
                      "Nested parallel with blocks are disallowed")

            self.state = 'parallel with'
            body = self.visitchild(node, 'body')
            self.state = None

            newnode.body = body
            return newnode
        elif self.parallel_directive:
            parallel_directive_class = self.get_directive_class_node(node)

            if not parallel_directive_class:
                # There was an error, stop here and now
                return None

            if parallel_directive_class is Nodes.ParallelWithBlockNode:
                error(node.pos, "The parallel directive must be called")
                return None

        self.visitchild(node, 'body')
        return node

    def visit_ForInStatNode(self, node):
        "Rewrite 'for i in cython.parallel.prange(...):'"
        self.visitchild(node, 'iterator')
        self.visitchild(node, 'target')

        in_prange = isinstance(node.iterator.sequence,
                               Nodes.ParallelRangeNode)
        previous_state = self.state

        if in_prange:
            # This will replace the entire ForInStatNode, so copy the
            # attributes
            parallel_range_node = node.iterator.sequence

            parallel_range_node.target = node.target
            parallel_range_node.body = node.body
            parallel_range_node.else_clause = node.else_clause

            node = parallel_range_node

            if not isinstance(node.target, ExprNodes.NameNode):
                error(node.target.pos,
                      "Can only iterate over an iteration variable")

            self.state = 'prange'

        self.visitchild(node, 'body')
        self.state = previous_state
        self.visitchild(node, 'else_clause')
        return node

    def visit(self, node):
        "Visit a node that may be None"
        if node is not None:
            return super().visit(node)


class WithTransform(VisitorTransform, SkipDeclarations):
    def visit_WithStatNode(self, node):
        self.visitchildren(node, ['body'])
        pos = node.pos
        is_async = node.is_async
        body, target, manager = node.body, node.target, node.manager
        manager = node.manager = ExprNodes.ProxyNode(manager)
        node.enter_call = ExprNodes.SimpleCallNode(
            pos, function=ExprNodes.AttributeNode(
                pos, obj=ExprNodes.CloneNode(manager),
                attribute=EncodedString('__aenter__' if is_async else '__enter__'),
                is_special_lookup=True),
            args=[],
            is_temp=True)

        if is_async:
            node.enter_call = ExprNodes.AwaitExprNode(pos, arg=node.enter_call)

        if target is not None:
            body = Nodes.StatListNode(
                pos, stats=[
                    Nodes.WithTargetAssignmentStatNode(
                        pos, lhs=target, with_node=node),
                    body])

        excinfo_target = ExprNodes.TupleNode(pos, slow=True, args=[
            ExprNodes.ExcValueNode(pos) for _ in range(3)])
        except_clause = Nodes.ExceptClauseNode(
            pos, body=Nodes.IfStatNode(
                pos, if_clauses=[
                    Nodes.IfClauseNode(
                        pos, condition=ExprNodes.NotNode(
                            pos, operand=ExprNodes.WithExitCallNode(
                                pos, with_stat=node,
                                test_if_run=False,
                                args=excinfo_target,
                                await_expr=ExprNodes.AwaitExprNode(pos, arg=None) if is_async else None)),
                        body=Nodes.ReraiseStatNode(pos),
                    ),
                ],
                else_clause=None),
            pattern=None,
            target=None,
            excinfo_target=excinfo_target,
        )

        node.body = Nodes.TryFinallyStatNode(
            pos, body=Nodes.TryExceptStatNode(
                pos, body=body,
                except_clauses=[except_clause],
                else_clause=None,
            ),
            finally_clause=Nodes.ExprStatNode(
                pos, expr=ExprNodes.WithExitCallNode(
                    pos, with_stat=node,
                    test_if_run=True,
                    args=ExprNodes.TupleNode(
                        pos, args=[ExprNodes.NoneNode(pos) for _ in range(3)]),
                    await_expr=ExprNodes.AwaitExprNode(pos, arg=None) if is_async else None)),
            handle_error_case=False,
        )
        return node

    def visit_ExprNode(self, node):
        # With statements are never inside expressions.
        return node

    visit_Node = VisitorTransform.recurse_to_children


class _GeneratorExpressionArgumentsMarker(TreeVisitor, SkipDeclarations):
    # called from "MarkClosureVisitor"
    def __init__(self, gen_expr):
        super().__init__()
        self.gen_expr = gen_expr

    def visit_ExprNode(self, node):
        if not node.is_literal:
            # Don't bother tagging literal nodes
            assert (not node.generator_arg_tag)  # nobody has tagged this first
            node.generator_arg_tag = self.gen_expr
        self.visitchildren(node)

    def visit_Node(self, node):
        # We're only interested in the expressions that make up the iterator sequence,
        # so don't go beyond ExprNodes (e.g. into ForFromStatNode).
        return

    def visit_GeneratorExpressionNode(self, node):
        node.generator_arg_tag = self.gen_expr
        # don't visit children, can't handle overlapping tags
        # (and assume generator expressions don't end up optimized out in a way
        #  that would require overlapping tags)


class _HandleGeneratorArguments(VisitorTransform, SkipDeclarations):
    # used from within CreateClosureClasses

    def __call__(self, node):
        from . import Visitor
        assert isinstance(node, ExprNodes.GeneratorExpressionNode)
        self.gen_node = node

        self.args = list(node.def_node.args)
        self.call_parameters = list(node.call_parameters)
        self.tag_count = 0
        self.substitutions = {}

        self.visitchildren(node)

        for k, v in self.substitutions.items():
            # doing another search for replacements here (at the end) allows us to sweep up
            # CloneNodes too (which are often generated by the optimizer)
            # (it could arguably be done more efficiently with a single traversal though)
            Visitor.recursively_replace_node(node, k, v)

        node.def_node.args = self.args
        node.call_parameters = self.call_parameters
        return node

    def visit_GeneratorExpressionNode(self, node):
        # a generator can also be substituted itself, so handle that case
        new_node = self._handle_ExprNode(node, do_visit_children=False)
        # However do not traverse into it. A new _HandleGeneratorArguments visitor will be used
        # elsewhere to do that.
        return node

    def _handle_ExprNode(self, node, do_visit_children):
        if (node.generator_arg_tag is not None and self.gen_node is not None and
                self.gen_node == node.generator_arg_tag):
            pos = node.pos
            # The reason for using ".x" as the name is that this is how CPython
            # tracks internal variables in loops (e.g.
            #  { locals() for v in range(10) }
            # will produce "v" and ".0"). We don't replicate this behaviour completely
            # but use it as a starting point
            name_source = self.tag_count
            self.tag_count += 1
            name = EncodedString(".{}".format(name_source))
            def_node = self.gen_node.def_node
            if not def_node.local_scope.lookup_here(name):
                from . import Symtab
                cname = EncodedString(Naming.genexpr_arg_prefix + Symtab.punycodify_name(str(name_source)))
                name_decl = Nodes.CNameDeclaratorNode(pos=pos, name=name)
                type = node.type

                # strip away cv types - they shouldn't be applied to the
                # function argument or to the closure struct.
                # It isn't obvious whether the right thing to do would be to capture by reference or by
                # value (C++ itself doesn't know either for lambda functions and forces a choice).
                # However, capture by reference involves converting to FakeReference which would require
                # re-analysing AttributeNodes. Therefore I've picked capture-by-value out of convenience
                # TODO - could probably be optimized by making the arg a reference but the closure not
                # (see https://github.com/cython/cython/issues/2468)
                type = PyrexTypes.remove_cv_ref(type, remove_fakeref=False)

                name_decl.type = type
                new_arg = Nodes.CArgDeclNode(pos=pos, declarator=name_decl,
                                                base_type=None, default=None, annotation=None)
                new_arg.name = name_decl.name
                new_arg.type = type

                self.args.append(new_arg)
                node.generator_arg_tag = None  # avoid the possibility of this being caught again
                self.call_parameters.append(node)
                new_arg.entry = def_node.declare_argument(def_node.local_scope, new_arg)
                new_arg.entry.cname = cname
                new_arg.entry.in_closure = True

            if do_visit_children:
                # now visit the Nodes's children (but remove self.gen_node to not to further
                # argument substitution)
                gen_node, self.gen_node = self.gen_node, None
                self.visitchildren(node)
                self.gen_node = gen_node

            # replace the node inside the generator with a looked-up name
            # (initialized_check can safely be False because the source variable will be checked
            # before it is captured if the check is required)
            name_node = ExprNodes.NameNode(pos, name=name, initialized_check=False)
            name_node.entry = self.gen_node.def_node.gbody.local_scope.lookup(name_node.name)
            name_node.type = name_node.entry.type
            self.substitutions[node] = name_node
            return name_node
        if do_visit_children:
            self.visitchildren(node)
        return node

    def visit_ExprNode(self, node):
        return self._handle_ExprNode(node, True)

    visit_Node = VisitorTransform.recurse_to_children


class DecoratorTransform(ScopeTrackingTransform, SkipDeclarations):
    """
    Transforms method decorators in cdef classes into nested calls or properties.

    Python-style decorator properties are transformed into a PropertyNode
    with up to the three getter, setter and deleter DefNodes.
    The functional style isn't supported yet.
    """
    _properties = None

    _map_property_attribute = {
        'getter': EncodedString('__get__'),
        'setter': EncodedString('__set__'),
        'deleter': EncodedString('__del__'),
    }.get

    def visit_CClassDefNode(self, node):
        if self._properties is None:
            self._properties = []
        self._properties.append({})
        node = super().visit_CClassDefNode(node)
        self._properties.pop()
        return node

    def visit_PropertyNode(self, node):
        # Low-level warning for other code until we can convert all our uses over.
        level = 2 if isinstance(node.pos[0], str) else 0
        warning(node.pos, "'property %s:' syntax is deprecated, use '@property'" % node.name, level)
        return node

    def visit_CFuncDefNode(self, node):
        node = self.visit_FuncDefNode(node)
        if not node.decorators:
            return node
        elif self.scope_type != 'cclass' or self.scope_node.visibility != "extern":
            # at the moment cdef functions are very restricted in what decorators they can take
            # so it's simple to test for the small number of allowed decorators....
            if not (len(node.decorators) == 1 and node.decorators[0].decorator.is_name and
                    node.decorators[0].decorator.name == "staticmethod"):
                error(node.decorators[0].pos, "Cdef functions cannot take arbitrary decorators.")
            return node

        ret_node = node
        decorator_node = self._find_property_decorator(node)
        if decorator_node:
            if decorator_node.decorator.is_name:
                name = node.declared_name()
                if name:
                    ret_node = self._add_property(node, name, decorator_node)
            else:
                error(decorator_node.pos, "C property decorator can only be @property")

        if node.decorators:
            return self._reject_decorated_property(node, node.decorators[0])
        return ret_node

    def visit_DefNode(self, node):
        scope_type = self.scope_type
        node = self.visit_FuncDefNode(node)
        if scope_type != 'cclass' or not node.decorators:
            return node

        # transform @property decorators
        decorator_node = self._find_property_decorator(node)
        if decorator_node is not None:
            decorator = decorator_node.decorator
            if decorator.is_name:
                return self._add_property(node, node.name, decorator_node)
            else:
                handler_name = self._map_property_attribute(decorator.attribute)
                if handler_name:
                    if decorator.obj.name != node.name:
                        # CPython does not generate an error or warning, but not something useful either.
                        error(decorator_node.pos,
                              "Mismatching property names, expected '%s', got '%s'" % (
                                  decorator.obj.name, node.name))
                    elif len(node.decorators) > 1:
                        return self._reject_decorated_property(node, decorator_node)
                    else:
                        return self._add_to_property(node, handler_name, decorator_node)

        # we clear node.decorators, so we need to set the
        # is_staticmethod/is_classmethod attributes now
        for decorator in node.decorators:
            func = decorator.decorator
            if func.is_name:
                node.is_classmethod |= func.name == 'classmethod'
                node.is_staticmethod |= func.name == 'staticmethod'

        # transform normal decorators
        decs = node.decorators
        node.decorators = None
        return self.chain_decorators(node, decs, node.name)

    def _find_property_decorator(self, node):
        properties = self._properties[-1]
        for decorator_node in node.decorators[::-1]:
            decorator = decorator_node.decorator
            if decorator.is_name and decorator.name == 'property':
                # @property
                return decorator_node
            elif decorator.is_attribute and decorator.obj.name in properties:
                # @prop.setter etc.
                return decorator_node
        return None

    @staticmethod
    def _reject_decorated_property(node, decorator_node):
        # restrict transformation to outermost decorator as wrapped properties will probably not work
        for deco in node.decorators:
            if deco != decorator_node:
                error(deco.pos, "Property methods with additional decorators are not supported")
        return node

    def _add_property(self, node, name, decorator_node):
        if len(node.decorators) > 1:
            return self._reject_decorated_property(node, decorator_node)
        node.decorators.remove(decorator_node)
        properties = self._properties[-1]
        is_cproperty = isinstance(node, Nodes.CFuncDefNode)
        body = Nodes.StatListNode(node.pos, stats=[node])
        if is_cproperty:
            if name in properties:
                error(node.pos, "C property redeclared")
            if 'inline' not in node.modifiers:
                error(node.pos, "C property method must be declared 'inline'")
            prop = Nodes.CPropertyNode(node.pos, doc=node.doc, name=name, body=body)
        elif name in properties:
            prop = properties[name]
            if prop.is_cproperty:
                error(node.pos, "C property redeclared")
            else:
                node.name = EncodedString("__get__")
                prop.pos = node.pos
                prop.doc = node.doc
                prop.body.stats = [node]
            return None
        else:
            node.name = EncodedString("__get__")
            prop = Nodes.PropertyNode(
                node.pos, name=name, doc=node.doc, body=body)
        properties[name] = prop
        return prop

    def _add_to_property(self, node, name, decorator):
        properties = self._properties[-1]
        prop = properties[node.name]
        if prop.is_cproperty:
            error(node.pos, "C property redeclared")
            return None
        node.name = name
        node.decorators.remove(decorator)
        stats = prop.body.stats
        for i, stat in enumerate(stats):
            if stat.name == name:
                stats[i] = node
                break
        else:
            stats.append(node)
        return None

    @staticmethod
    def chain_decorators(node, decorators, name):
        """
        Decorators are applied directly in DefNode and PyClassDefNode to avoid
        reassignments to the function/class name - except for cdef class methods.
        For those, the reassignment is required as methods are originally
        defined in the PyMethodDef struct.

        The IndirectionNode allows DefNode to override the decorator.
        """
        decorator_result = ExprNodes.NameNode(node.pos, name=name)
        for decorator in decorators[::-1]:
            decorator_result = ExprNodes.SimpleCallNode(
                decorator.pos,
                function=decorator.decorator,
                args=[decorator_result])

        name_node = ExprNodes.NameNode(node.pos, name=name)
        reassignment = Nodes.SingleAssignmentNode(
            node.pos,
            lhs=name_node,
            rhs=decorator_result)

        reassignment = Nodes.IndirectionNode([reassignment])
        node.decorator_indirection = reassignment
        return [node, reassignment]


class CnameDirectivesTransform(CythonTransform, SkipDeclarations):
    """
    Only part of the CythonUtilityCode pipeline. Must be run before
    DecoratorTransform in case this is a decorator for a cdef class.
    It filters out @cname('my_cname') decorators and rewrites them to
    CnameDecoratorNodes.
    """

    def handle_function(self, node):
        if not getattr(node, 'decorators', None):
            return self.visit_Node(node)

        for i, decorator in enumerate(node.decorators):
            decorator = decorator.decorator

            if (isinstance(decorator, ExprNodes.CallNode) and
                    decorator.function.is_name and
                    decorator.function.name == 'cname'):
                args, kwargs = decorator.explicit_args_kwds()

                if kwargs:
                    raise AssertionError(
                            "cname decorator does not take keyword arguments")

                if len(args) != 1:
                    raise AssertionError(
                            "cname decorator takes exactly one argument")

                if not (args[0].is_literal and args[0].type is Builtin.unicode_type):
                    raise AssertionError(
                            "argument to cname decorator must be a string literal")

                cname = args[0].compile_time_value(None)
                del node.decorators[i]
                node = Nodes.CnameDecoratorNode(pos=node.pos, node=node,
                                                cname=cname)
                break

        return self.visit_Node(node)

    visit_FuncDefNode = handle_function
    visit_CClassDefNode = handle_function
    visit_CEnumDefNode = handle_function
    visit_CStructOrUnionDefNode = handle_function
    visit_CVarDefNode = handle_function


class ForwardDeclareTypes(CythonTransform):
    """
    Declare all global cdef names that we allow referencing in other places,
    before declaring everything (else) in source code order.
    """

    def visit_CompilerDirectivesNode(self, node):
        env = self.module_scope
        old = env.directives
        env.directives = node.directives
        self.visitchildren(node)
        env.directives = old
        return node

    def visit_ModuleNode(self, node):
        self.module_scope = node.scope
        self.module_scope.directives = node.directives
        self.visitchildren(node)
        return node

    def visit_CDefExternNode(self, node):
        old_cinclude_flag = self.module_scope.in_cinclude
        self.module_scope.in_cinclude = 1
        self.visitchildren(node)
        self.module_scope.in_cinclude = old_cinclude_flag
        return node

    def visit_CEnumDefNode(self, node):
        node.declare(self.module_scope)
        return node

    def visit_CStructOrUnionDefNode(self, node):
        if node.name not in self.module_scope.entries:
            node.declare(self.module_scope)
        return node

    def visit_CClassDefNode(self, node):
        if node.class_name not in self.module_scope.entries:
            node.declare(self.module_scope)
        # Expand fused methods of .pxd declared types to construct the final vtable order.
        type = self.module_scope.entries[node.class_name].type
        if type is not None and type.is_extension_type and not type.is_builtin_type and type.scope:
            scope = type.scope
            for entry in scope.cfunc_entries:
                if entry.type and entry.type.is_fused:
                    entry.type.get_all_specialized_function_types()
        return node

    def visit_FuncDefNode(self, node):
        # no traversal needed
        return node

    def visit_PyClassDefNode(self, node):
        # no traversal needed
        return node


class AnalyseDeclarationsTransform(EnvTransform):

    basic_property = TreeFragment("""
property NAME:
    def __get__(self):
        return ATTR
    def __set__(self, value):
        ATTR = value
    """, level='c_class', pipeline=[NormalizeTree(None)])
    basic_pyobject_property = TreeFragment("""
property NAME:
    def __get__(self):
        return ATTR
    def __set__(self, value):
        ATTR = value
    def __del__(self):
        ATTR = None
    """, level='c_class', pipeline=[NormalizeTree(None)])
    basic_property_ro = TreeFragment("""
property NAME:
    def __get__(self):
        return ATTR
    """, level='c_class', pipeline=[NormalizeTree(None)])

    struct_or_union_wrapper = TreeFragment("""
cdef class NAME:
    cdef TYPE value
    def __init__(self, MEMBER=None):
        cdef int count
        count = 0
        INIT_ASSIGNMENTS
        if IS_UNION and count > 1:
            raise ValueError, "At most one union member should be specified."
    def __str__(self):
        return STR_FORMAT % MEMBER_TUPLE
    def __repr__(self):
        return REPR_FORMAT % MEMBER_TUPLE
    """, pipeline=[NormalizeTree(None)])

    init_assignment = TreeFragment("""
if VALUE is not None:
    ATTR = VALUE
    count += 1
    """, pipeline=[NormalizeTree(None)])

    fused_function = None
    in_lambda = 0

    def __call__(self, root):
        # needed to determine if a cdef var is declared after it's used.
        self.seen_vars_stack = []
        self.fused_error_funcs = set()
        super_class = super()
        self._super_visit_FuncDefNode = super_class.visit_FuncDefNode
        return super_class.__call__(root)

    def visit_NameNode(self, node):
        self.seen_vars_stack[-1].add(node.name)
        return node

    def visit_ModuleNode(self, node):
        # Pickling support requires injecting module-level nodes.
        self.extra_module_declarations = []
        self.seen_vars_stack.append(set())
        node.analyse_declarations(self.current_env())
        self.visitchildren(node)
        self.seen_vars_stack.pop()
        node.body.stats.extend(self.extra_module_declarations)
        return node

    def visit_LambdaNode(self, node):
        self.in_lambda += 1
        node.analyse_declarations(self.current_env())
        self.visitchildren(node)
        self.in_lambda -= 1
        return node

    def visit_CClassDefNode(self, node):
        node = self.visit_ClassDefNode(node)
        if node.scope and 'dataclasses.dataclass' in node.scope.directives:
            from .Dataclass import handle_cclass_dataclass
            handle_cclass_dataclass(node, node.scope.directives['dataclasses.dataclass'], self)
        if node.scope and node.scope.implemented and node.body:
            stats = []
            for entry in node.scope.var_entries:
                if entry.needs_property:
                    property = self.create_Property(entry)
                    property.analyse_declarations(node.scope)
                    self.visit(property)
                    stats.append(property)
            if stats:
                node.body.stats += stats
            if (node.visibility != 'extern'
                    and not node.scope.lookup('__reduce__')
                    and not node.scope.lookup('__reduce_ex__')):
                self._inject_pickle_methods(node)
        return node

    def _inject_pickle_methods(self, node):
        env = self.current_env()
        if node.scope.directives['auto_pickle'] is False:   # None means attempt it.
            # Old behavior of not doing anything.
            return
        auto_pickle_forced = node.scope.directives['auto_pickle'] is True

        all_members = []
        cls = node.entry.type
        cinit = None
        inherited_reduce = None
        while cls is not None:
            all_members.extend(e for e in cls.scope.var_entries if e.name not in ('__weakref__', '__dict__'))
            cinit = cinit or cls.scope.lookup('__cinit__')
            inherited_reduce = inherited_reduce or cls.scope.lookup('__reduce__') or cls.scope.lookup('__reduce_ex__')
            cls = cls.base_type
        all_members.sort(key=lambda e: e.name)

        if inherited_reduce:
            # This is not failsafe, as we may not know whether a cimported class defines a __reduce__.
            # This is why we define __reduce_cython__ and only replace __reduce__
            # (via ExtensionTypes.SetupReduce utility code) at runtime on class creation.
            return

        non_py = [
            e for e in all_members
            if not e.type.is_pyobject and (not e.type.can_coerce_to_pyobject(env)
                                           or not e.type.can_coerce_from_pyobject(env))
        ]

        structs = [e for e in all_members if e.type.is_struct_or_union]

        if cinit or non_py or (structs and not auto_pickle_forced):
            if cinit:
                # TODO(robertwb): We could allow this if __cinit__ has no require arguments.
                msg = 'no default __reduce__ due to non-trivial __cinit__'
            elif non_py:
                msg = "%s cannot be converted to a Python object for pickling" % ','.join("self.%s" % e.name for e in non_py)
            else:
                # Extern structs may be only partially defined.
                # TODO(robertwb): Limit the restriction to extern
                # (and recursively extern-containing) structs.
                msg = ("Pickling of struct members such as %s must be explicitly requested "
                       "with @auto_pickle(True)" % ','.join("self.%s" % e.name for e in structs))

            if auto_pickle_forced:
                error(node.pos, msg)

            pickle_func = TreeFragment("""
                def __reduce_cython__(self):
                    raise TypeError, "%(msg)s"
                def __setstate_cython__(self, __pyx_state):
                    raise TypeError, "%(msg)s"
                """ % {'msg': msg},
                level='c_class', pipeline=[NormalizeTree(None)]).substitute({})
            pickle_func.analyse_declarations(node.scope)
            self.visit(pickle_func)
            node.body.stats.append(pickle_func)

        else:
            for e in all_members:
                if not e.type.is_pyobject:
                    e.type.create_to_py_utility_code(env)
                    e.type.create_from_py_utility_code(env)

            all_members_names = [e.name for e in all_members]
            assignments = '; '.join([
                '__pyx_result.%s = __pyx_state[%s]' % (v, ix)
                for ix, v in enumerate(all_members_names)
            ])
            checksums = _calculate_pickle_checksums(all_members_names)
            if len(checksums) != 3:
                # If we don't have enough checksums to call the check function, we just repeat the last one.
                checksums = (checksums + [checksums[-1] * 2])[:3]

            unpickle_func_name = f'__pyx_unpickle_{node.punycode_class_name}'
            num_members = len(all_members_names)

            env.use_utility_code(Code.UtilityCode.load_cached("UpdateUnpickledDict", "ExtensionTypes.c"))

            # TODO(robertwb): Move the state into the third argument
            # so it can be pickled *after* self is memoized.
            unpickle_code = f"""
                cdef extern from *:
                    int __Pyx_CheckUnpickleChecksum(long, long, long, long, const char*) except -1
                    int __Pyx_UpdateUnpickledDict(object, object, Py_ssize_t) except -1

                def {unpickle_func_name}(__pyx_type, long __pyx_checksum, tuple __pyx_state):
                    cdef object __pyx_result
                    __Pyx_CheckUnpickleChecksum(__pyx_checksum, {', '.join(checksums)}, {', '.join(all_members_names).encode('UTF-8')!r})
                    __pyx_result = {node.class_name}.__new__(__pyx_type)
                    if __pyx_state is not None:
                        {unpickle_func_name}__set_state(<{node.class_name}> __pyx_result, __pyx_state)
                    return __pyx_result

                cdef {unpickle_func_name}__set_state({node.class_name} __pyx_result, __pyx_state: tuple):
                    {assignments}
                    __Pyx_UpdateUnpickledDict(__pyx_result, __pyx_state, {num_members:d})
            """

            env.use_utility_code(Code.UtilityCode.load_cached("CheckUnpickleChecksum", "ExtensionTypes.c"))

            unpickle_func = TreeFragment(unpickle_code, level='module', pipeline=[NormalizeTree(None)]).substitute({})
            unpickle_func.analyse_declarations(node.entry.scope)

            self.visit(unpickle_func)
            self.extra_module_declarations.append(unpickle_func)

            members = ', '.join(f'self.{v}' for v in all_members_names) + (',' if len(all_members_names) == 1 else '')
            # Even better, we could check PyType_IS_GC.
            any_notnone_members = ' or '.join([f'self.{e.name} is not None' for e in all_members if e.type.is_pyobject] or ['False'])

            pickle_code = f"""
                def __reduce_cython__(self):
                    cdef tuple state
                    cdef object _dict
                    cdef bint use_setstate
                    state = ({members})
                    _dict = getattr(self, '__dict__', None)
                    if _dict is not None and _dict:
                        state += (_dict,)
                        use_setstate = True
                    else:
                        use_setstate = {any_notnone_members}
                    if use_setstate:
                        return {unpickle_func_name}, (type(self), {checksums[0]}, None), state
                    else:
                        return {unpickle_func_name}, (type(self), {checksums[0]}, state)

                def __setstate_cython__(self, __pyx_state):
                    {unpickle_func_name}__set_state(self, __pyx_state)
            """

            pickle_func = TreeFragment(pickle_code, level='c_class', pipeline=[NormalizeTree(None)]).substitute({})
            pickle_func.analyse_declarations(node.scope)

            self.enter_scope(node, node.scope)  # functions should be visited in the class scope
            self.visit(pickle_func)
            self.exit_scope()
            node.body.stats.append(pickle_func)

    def _handle_fused_def_decorators(self, old_decorators, env, node):
        """
        Create function calls to the decorators and reassignments to
        the function.
        """
        # Delete staticmethod and classmethod decorators, this is
        # handled directly by the fused function object.
        decorators = []
        for decorator in old_decorators:
            func = decorator.decorator
            if (not func.is_name or
                    func.name not in ('staticmethod', 'classmethod') or
                    env.lookup_here(func.name)):
                # not a static or classmethod
                decorators.append(decorator)

        if decorators:
            transform = DecoratorTransform(self.context)
            def_node = node.node
            _, reassignments = transform.chain_decorators(
                def_node, decorators, def_node.name)
            reassignments.analyse_declarations(env)
            node = [node, reassignments]

        return node

    def _handle_def(self, decorators, env, node):
        "Handle def or cpdef fused functions"
        # Create PyCFunction nodes for each specialization
        node.stats.insert(0, node.py_func)
        self.visitchild(node, 'py_func')
        node.update_fused_defnode_entry(env)
        # For the moment, fused functions do not support METH_FASTCALL
        node.py_func.entry.signature.use_fastcall = False
        pycfunc = ExprNodes.PyCFunctionNode.from_defnode(node.py_func, binding=True)
        pycfunc = ExprNodes.ProxyNode(pycfunc.coerce_to_temp(env))
        node.resulting_fused_function = pycfunc
        # Create assignment node for our def function
        node.fused_func_assignment = self._create_assignment(
            node.py_func, ExprNodes.CloneNode(pycfunc), env)

        if decorators:
            node = self._handle_fused_def_decorators(decorators, env, node)

        return node

    def _create_fused_function(self, env, node):
        "Create a fused function for a DefNode with fused arguments"
        from . import FusedNode

        if self.fused_function or self.in_lambda:
            if self.fused_function not in self.fused_error_funcs:
                if self.in_lambda:
                    error(node.pos, "Fused lambdas not allowed")
                else:
                    error(node.pos, "Cannot nest fused functions")

            self.fused_error_funcs.add(self.fused_function)

            node.body = Nodes.PassStatNode(node.pos)
            for arg in node.args:
                if arg.type.is_fused:
                    arg.type = arg.type.get_fused_types()[0]

            return node

        decorators = getattr(node, 'decorators', None)
        node = FusedNode.FusedCFuncDefNode(node, env)
        self.fused_function = node
        self.visitchildren(node)
        self.fused_function = None
        if node.py_func:
            node = self._handle_def(decorators, env, node)

        return node

    def _handle_fused(self, node):
        if node.is_generator and node.has_fused_arguments:
            error(node.pos, "Fused generators not supported")
            node.has_fused_arguments = False
            node.gbody.body = Nodes.StatListNode(node.pos, stats=[])

        return node.has_fused_arguments

    def visit_FuncDefNode(self, node):
        """
        Analyse a function and its body, as that hasn't happened yet.  Also
        analyse the directive_locals set by @cython.locals().

        Then, if we are a function with fused arguments, replace the function
        (after it has declared itself in the symbol table!) with a
        FusedCFuncDefNode, and analyse its children (which are in turn normal
        functions). If we're a normal function, just analyse the body of the
        function.
        """
        env = self.current_env()

        self.seen_vars_stack.append(set())
        lenv = node.local_scope
        node.declare_arguments(lenv)

        # @cython.locals(...)
        for var, type_node in node.directive_locals.items():
            if not lenv.lookup_here(var):   # don't redeclare args
                type = type_node.analyse_as_type(lenv)
                if type and type.is_fused and lenv.fused_to_specific:
                    type = type.specialize(lenv.fused_to_specific)
                if type:
                    lenv.declare_var(var, type, type_node.pos)
                else:
                    error(type_node.pos, "Not a type")

        if self._handle_fused(node):
            node = self._create_fused_function(env, node)
        else:
            node.body.analyse_declarations(lenv)
            node = self._super_visit_FuncDefNode(node)

        self.seen_vars_stack.pop()

        if "ufunc" in lenv.directives:
            from . import UFuncs
            return UFuncs.convert_to_ufunc(node)
        return node

    def visit_DefNode(self, node):
        node = self.visit_FuncDefNode(node)
        if not isinstance(node, Nodes.DefNode):
            return node
        env = self.current_env()
        if node.code_object is None:
            node.code_object = ExprNodes.CodeObjectNode(node)
            node.code_object.analyse_declarations(env)
        if node.fused_py_func or node.is_generator_body:
            return node
        if not node.needs_assignment_synthesis(env):
            return node
        return [node, self._synthesize_assignment(node, env)]

    def visit_CFuncDefNode(self, node):
        if node.code_object is None and node.py_func is None:
            node.code_object = ExprNodes.CodeObjectNode.for_cfunc(node)
            node.code_object.analyse_declarations(self.current_env())
        return self.visit_FuncDefNode(node)

    def visit_GeneratorBodyDefNode(self, node):
        return self.visit_FuncDefNode(node)

    def visit_GeneratorDefNode(self, node):
        # The generator body should use the same code object as the (user facing) generator function that creates it.
        result = self.visit_DefNode(node)
        # 'result' will usually be a list of statements, but we still have the original node.
        node.gbody.code_object = node.code_object
        return result

    def _synthesize_assignment(self, node, env):
        # Synthesize assignment node and put it right after defnode
        genv = env
        while genv.is_py_class_scope or genv.is_c_class_scope:
            genv = genv.outer_scope

        binding = env.is_py_class_scope or self.current_directives.get('binding')
        if genv.is_closure_scope:
            rhs = node.py_cfunc_node = ExprNodes.InnerFunctionNode.from_defnode(node, binding)
        else:
            rhs = ExprNodes.PyCFunctionNode.from_defnode(node, binding)

        node.is_cyfunction = rhs.binding
        return self._create_assignment(node, rhs, env)

    def _create_assignment(self, def_node, rhs, env):
        if def_node.decorators:
            for decorator in def_node.decorators[::-1]:
                rhs = ExprNodes.SimpleCallNode(
                    decorator.pos,
                    function = decorator.decorator,
                    args = [rhs])
            def_node.decorators = None

        assmt = Nodes.SingleAssignmentNode(
            def_node.pos,
            lhs=ExprNodes.NameNode(def_node.pos, name=def_node.name),
            rhs=rhs)
        assmt.analyse_declarations(env)
        return assmt

    def visit_func_outer_attrs(self, node):
        # any names in the outer attrs should not be looked up in the function "seen_vars_stack"
        stack = self.seen_vars_stack.pop()
        super().visit_func_outer_attrs(node)
        self.seen_vars_stack.append(stack)

    def visit_ScopedExprNode(self, node):
        env = self.current_env()
        node.analyse_declarations(env)
        # the node may or may not have a local scope
        if node.expr_scope:
            self.seen_vars_stack.append(set(self.seen_vars_stack[-1]))
            self.enter_scope(node, node.expr_scope)
            node.analyse_scoped_declarations(node.expr_scope)
            self.visitchildren(node)
            self.exit_scope()
            self.seen_vars_stack.pop()
        else:

            node.analyse_scoped_declarations(env)
            self.visitchildren(node)
        return node

    def visit_TempResultFromStatNode(self, node):
        self.visitchildren(node)
        node.analyse_declarations(self.current_env())
        return node

    def visit_CppClassNode(self, node):
        if node.visibility == 'extern':
            return None
        else:
            return self.visit_ClassDefNode(node)

    def visit_CStructOrUnionDefNode(self, node):
        # Create a wrapper node if needed.
        # We want to use the struct type information (so it can't happen
        # before this phase) but also create new objects to be declared
        # (so it can't happen later).
        # Note that we don't return the original node, as it is
        # never used after this phase.
        if True:  # private (default)
            return None

        self_value = ExprNodes.AttributeNode(
            pos = node.pos,
            obj = ExprNodes.NameNode(pos=node.pos, name="self"),
            attribute = EncodedString("value"))
        var_entries = node.entry.type.scope.var_entries
        attributes = []
        for entry in var_entries:
            attributes.append(ExprNodes.AttributeNode(pos = entry.pos,
                                                      obj = self_value,
                                                      attribute = entry.name))
        # __init__ assignments
        init_assignments = []
        for entry, attr in zip(var_entries, attributes):
            # TODO: branch on visibility
            init_assignments.append(
                self.init_assignment.substitute(
                    {
                        "VALUE": ExprNodes.NameNode(entry.pos, name = entry.name),
                        "ATTR": attr,
                    },
                    pos=entry.pos,
                )
            )

        # create the class
        str_format = "%s(%s)" % (node.entry.type.name, ("%s, " * len(attributes))[:-2])
        wrapper_class = self.struct_or_union_wrapper.substitute({
            "INIT_ASSIGNMENTS": Nodes.StatListNode(node.pos, stats = init_assignments),
            "IS_UNION": ExprNodes.BoolNode(node.pos, value = not node.entry.type.is_struct),
            "MEMBER_TUPLE": ExprNodes.TupleNode(node.pos, args=attributes),
            "STR_FORMAT": ExprNodes.UnicodeNode(node.pos, value = EncodedString(str_format)),
            "REPR_FORMAT": ExprNodes.UnicodeNode(node.pos, value = EncodedString(str_format.replace("%s", "%r"))),
        }, pos = node.pos).stats[0]
        wrapper_class.class_name = node.name
        wrapper_class.shadow = True
        class_body = wrapper_class.body.stats

        # fix value type
        assert isinstance(class_body[0].base_type, Nodes.CSimpleBaseTypeNode)
        class_body[0].base_type.name = node.name

        # fix __init__ arguments
        init_method = class_body[1]
        assert isinstance(init_method, Nodes.DefNode) and init_method.name == '__init__'
        arg_template = init_method.args[1]
        if not node.entry.type.is_struct:
            arg_template.kw_only = True
        del init_method.args[1]
        for entry, attr in zip(var_entries, attributes):
            arg = copy.deepcopy(arg_template)
            arg.declarator.name = entry.name
            init_method.args.append(arg)

        # setters/getters
        for entry, attr in zip(var_entries, attributes):
            # TODO: branch on visibility
            if entry.type.is_pyobject:
                template = self.basic_pyobject_property
            else:
                template = self.basic_property
            property = template.substitute(
                {
                    "ATTR": attr,
                },
                pos=entry.pos,
            ).stats[0]
            property.name = entry.name
            wrapper_class.body.stats.append(property)

        wrapper_class.analyse_declarations(self.current_env())
        return self.visit_CClassDefNode(wrapper_class)

    # Some nodes are no longer needed after declaration
    # analysis and can be dropped. The analysis was performed
    # on these nodes in a separate recursive process from the
    # enclosing function or module, so we can simply drop them.
    def visit_CDeclaratorNode(self, node):
        # necessary to ensure that all CNameDeclaratorNodes are visited.
        self.visitchildren(node)
        return node

    def visit_CTypeDefNode(self, node):
        return node

    def visit_CBaseTypeNode(self, node):
        return None

    def visit_CEnumDefNode(self, node):
        if node.visibility == 'public':
            return node
        else:
            return None

    def visit_CNameDeclaratorNode(self, node):
        if node.name in self.seen_vars_stack[-1]:
            entry = self.current_env().lookup(node.name)
            if (entry is None or entry.visibility != 'extern'
                    and not entry.scope.is_c_class_scope):
                error(node.pos, "cdef variable '%s' declared after it is used" % node.name)
        self.visitchildren(node)
        return node

    def visit_CVarDefNode(self, node):
        # to ensure all CNameDeclaratorNodes are visited.
        self.visitchildren(node)
        return None

    def visit_CnameDecoratorNode(self, node):
        child_node = self.visitchild(node, 'node')
        if not child_node:
            return None
        if type(child_node) is list:  # Assignment synthesized
            node.node = child_node[0]
            return [node] + child_node[1:]
        return node

    def create_Property(self, entry):
        if entry.visibility == 'public':
            if entry.type.is_pyobject:
                template = self.basic_pyobject_property
            else:
                template = self.basic_property
        elif entry.visibility == 'readonly':
            template = self.basic_property_ro
        property = template.substitute(
            {
                "ATTR": ExprNodes.AttributeNode(pos=entry.pos,
                                                obj=ExprNodes.NameNode(pos=entry.pos, name="self"),
                                                attribute=entry.name),
            },
            pos=entry.pos,
        ).stats[0]
        property.name = entry.name
        property.doc = entry.doc
        return property

    def visit_AssignmentExpressionNode(self, node):
        self.visitchildren(node)
        node.analyse_declarations(self.current_env())
        return node


def _calculate_pickle_checksums(member_names):
    # Cython 0.x used MD5 for the checksum, which a few Python installations remove for security reasons.
    # SHA-256 should be ok for years to come, but early Cython 3.0 alpha releases used SHA-1,
    # which may not be.
    member_names_string = ' '.join(member_names).encode('utf-8')
    hash_kwargs = {'usedforsecurity': False} if sys.version_info >= (3, 9) else {}
    checksums = []
    for algo_name in ['sha256', 'sha1', 'md5']:
        try:
            mkchecksum = getattr(hashlib, algo_name)
            checksum = mkchecksum(member_names_string, **hash_kwargs).hexdigest()
        except (AttributeError, ValueError):
            # The algorithm (i.e. MD5) might not be there at all, or might be blocked at runtime.
            continue
        checksums.append('0x' + checksum[:7])
    return checksums


class CalculateQualifiedNamesTransform(EnvTransform):
    """
    Calculate and store the '__qualname__' and the global
    module name on some nodes.
    """
    needs_qualname_assignment = False
    needs_module_assignment = False

    def visit_ModuleNode(self, node):
        self.module_name = self.global_scope().qualified_name
        self.qualified_name = []
        _super = super()
        self._super_visit_FuncDefNode = _super.visit_FuncDefNode
        self._super_visit_ClassDefNode = _super.visit_ClassDefNode
        self.visitchildren(node)
        return node

    def _set_qualname(self, node, name=None):
        if name:
            qualname = self.qualified_name[:]
            qualname.append(name)
        else:
            qualname = self.qualified_name
        node.qualname = EncodedString('.'.join(qualname))
        node.module_name = self.module_name

    def _append_entry(self, entry):
        if entry.is_pyglobal and not entry.is_pyclass_attr:
            self.qualified_name = [entry.name]
        else:
            self.qualified_name.append(entry.name)

    def visit_ClassNode(self, node):
        self._set_qualname(node, node.name)
        self.visitchildren(node)
        return node

    def visit_PyClassNamespaceNode(self, node):
        # class name was already added by parent node
        self._set_qualname(node)
        self.visitchildren(node)
        return node

    def visit_PyCFunctionNode(self, node):
        orig_qualified_name = self.qualified_name[:]
        if node.def_node.is_wrapper and self.qualified_name and self.qualified_name[-1] == '<locals>':
            self.qualified_name.pop()
            self._set_qualname(node)
        else:
            self._set_qualname(node, node.def_node.name)
        self.visitchildren(node)
        self.qualified_name = orig_qualified_name
        return node

    def visit_DefNode(self, node):
        if node.is_wrapper and self.qualified_name:
            assert self.qualified_name[-1] == '<locals>', self.qualified_name
            orig_qualified_name = self.qualified_name[:]
            self.qualified_name.pop()
            self._set_qualname(node)
            self._super_visit_FuncDefNode(node)
            self.qualified_name = orig_qualified_name
        else:
            self._set_qualname(node, node.name)
            self.visit_FuncDefNode(node)
        return node

    def visit_FuncDefNode(self, node):
        orig_qualified_name = self.qualified_name[:]
        if getattr(node, 'name', None) == '<lambda>':
            self.qualified_name.append('<lambda>')
        else:
            self._append_entry(node.entry)
        self.qualified_name.append('<locals>')
        self._super_visit_FuncDefNode(node)
        self.qualified_name = orig_qualified_name
        return node

    def generate_assignment(self, node, name, value):
        entry = node.scope.lookup_here(name)
        lhs = ExprNodes.NameNode(
            node.pos,
            name=EncodedString(name),
            entry=entry,
            is_target=True)
        rhs = ExprNodes.UnicodeNode(node.pos, value=value)
        node.body.stats.insert(0, Nodes.SingleAssignmentNode(
            node.pos,
            lhs=lhs,
            rhs=rhs,
        ).analyse_expressions(self.current_env()))

    def visit_ClassDefNode(self, node):
        orig_needs_qualname_assignment = self.needs_qualname_assignment
        self.needs_qualname_assignment = False
        orig_needs_module_assignment = self.needs_module_assignment
        self.needs_module_assignment = False
        orig_qualified_name = self.qualified_name[:]
        entry = (getattr(node, 'entry', None) or             # PyClass
                 self.current_env().lookup_here(node.target.name))  # CClass
        self._append_entry(entry)
        self._super_visit_ClassDefNode(node)
        if self.needs_qualname_assignment:
            self.generate_assignment(node, "__qualname__",
                                     EncodedString(".".join(self.qualified_name)))
        if self.needs_module_assignment:
            self.generate_assignment(node, "__module__",
                                     EncodedString(self.module_name))
        self.qualified_name = orig_qualified_name
        self.needs_qualname_assignment = orig_needs_qualname_assignment
        self.needs_module_assignment = orig_needs_module_assignment
        return node

    def visit_NameNode(self, node):
        scope = self.current_env()
        if scope.is_c_class_scope:
            # unlike for a PyClass scope, these attributes aren't defined in the
            # dictionary when the class definition is executed, therefore we ask
            # the compiler to generate an assignment to them at the start of the
            # body.
            # NOTE: this doesn't put them in locals()
            if node.name == "__qualname__":
                self.needs_qualname_assignment = True
            elif node.name == "__module__":
                self.needs_module_assignment = True
        return node


class AnalyseExpressionsTransform(CythonTransform):

    def visit_ModuleNode(self, node):
        node.scope.infer_types()
        node.body = node.body.analyse_expressions(node.scope)
        self.positions = [{node.pos}]
        self.visitchildren(node)
        self._build_positions(node)
        return node

    def visit_FuncDefNode(self, node):
        node.local_scope.infer_types()
        node.body = node.body.analyse_expressions(node.local_scope)
        self.positions[-1].add(node.pos)

        if node.is_wrapper:
            # Share positions between function and Python wrapper.
            local_positions = self.positions[-1]
        else:
            local_positions = {node.pos}
        self.positions.append(local_positions)

        self.visitchildren(node)
        self._build_positions(node)
        return node

    def visit_ScopedExprNode(self, node):
        if node.has_local_scope:
            node.expr_scope.infer_types()
            node = node.analyse_scoped_expressions(node.expr_scope)
        self.visit_ExprNode(node)
        return node

    def visit_IndexNode(self, node):
        """
        Replace index nodes used to specialize cdef functions with fused
        argument types with the Attribute- or NameNode referring to the
        function. We then need to copy over the specialization properties to
        the attribute or name node.

        Because the indexing might be a Python indexing operation on a fused
        function, or (usually) a Cython indexing operation, we need to
        re-analyse the types.
        """
        self.visit_ExprNode(node)
        if node.is_fused_index and not node.type.is_error:
            node = node.base
        return node

    # Build the line table according to PEP-626.
    # We mostly just do this here to avoid yet another transform traversal.

    def visit_ExprNode(self, node):
        self.positions[-1].add(node.pos)
        self.visitchildren(node)
        return node

    def visit_StatNode(self, node):
        self.positions[-1].add(node.pos)
        self.visitchildren(node)
        return node

    def _build_positions(self, func_node):
        """
        Build the PEP-626 line table and "bytecode-to-position" mapping used for CodeObjects.
        """
        # Code can originate from different source files and string code fragments, even within a single function.
        # Thus, it's not completely correct to just ignore the source files when sorting the line numbers,
        # but it also doesn't hurt much for the moment. Eventually, we might need different CodeObjects
        # even within a single function if it uses code from different sources / line number ranges.
        positions: list = sorted(
            self.positions.pop(),
            key=itemgetter(1, 2),  # (line, column)
            # Build ranges backwards to know the end column before we see the start column in the same line.
            reverse=True,
        )

        next_line = -1
        next_column_in_line = 0

        ranges = []
        for _, line, start_column in positions:
            ranges.append((line, line, start_column, next_column_in_line if line == next_line else start_column + 1))
            next_line, next_column_in_line = line, start_column

        ranges.reverse()
        func_node.node_positions = ranges

        positions.reverse()
        i: cython.Py_ssize_t
        func_node.local_scope.node_positions_to_offset = {
            position: i
            for i, position in enumerate(positions)
        }


class FindInvalidUseOfFusedTypes(TreeVisitor):

    def __call__(self, tree):
        self._in_fused_function = False
        self.visit(tree)
        return tree

    def visit_Node(self, node):
        self.visitchildren(node)

    def visit_FuncDefNode(self, node):
        outer_status = self._in_fused_function
        self._in_fused_function = node.has_fused_arguments

        if not self._in_fused_function:
            # Errors related to use in functions with fused args will already
            # have been detected.
            if not node.is_generator_body and node.return_type.is_fused:
                error(node.pos, "Return type is not specified as argument type")

        self.visitchildren(node)
        self._in_fused_function = outer_status

    def visit_ExprNode(self, node):
        if not self._in_fused_function and node.type and node.type.is_fused:
            error(node.pos, "Invalid use of fused types, type cannot be specialized")
            # Errors in subtrees are likely related, so do not recurse.
        else:
            self.visitchildren(node)


class ExpandInplaceOperators(EnvTransform):

    def visit_InPlaceAssignmentNode(self, node):
        lhs = node.lhs
        rhs = node.rhs
        if lhs.type.is_cpp_class:
            # No getting around this exact operator here.
            return node
        if isinstance(lhs, ExprNodes.BufferIndexNode):
            # There is code to handle this case in InPlaceAssignmentNode
            return node

        env = self.current_env()
        def side_effect_free_reference(node, setting=False):
            if node.is_name:
                return node, []
            elif node.type.is_pyobject and not setting:
                node = LetRefNode(node)
                return node, [node]
            elif node.is_subscript:
                base, temps = side_effect_free_reference(node.base)
                index = LetRefNode(node.index)
                return ExprNodes.IndexNode(node.pos, base=base, index=index), temps + [index]
            elif node.is_attribute:
                obj, temps = side_effect_free_reference(node.obj, setting=setting)
                return ExprNodes.AttributeNode(node.pos, obj=obj, attribute=node.attribute), temps
            elif isinstance(node, ExprNodes.BufferIndexNode):
                raise ValueError("Don't allow things like attributes of buffer indexing operations")
            else:
                node = LetRefNode(node)
                return node, [node]
        try:
            lhs, let_ref_nodes = side_effect_free_reference(lhs, setting=True)
        except ValueError:
            return node
        dup = lhs.__class__(**lhs.__dict__)
        binop = ExprNodes.binop_node(node.pos,
                                     operator = node.operator,
                                     operand1 = dup,
                                     operand2 = rhs,
                                     inplace=True)
        # Manually analyse types for new node.
        lhs.is_target = True
        lhs = lhs.analyse_target_types(env)
        dup.analyse_types(env)  # FIXME: no need to reanalyse the copy, right?
        binop.analyse_operation(env)
        node = Nodes.SingleAssignmentNode(
            node.pos,
            lhs = lhs,
            rhs=binop.coerce_to(lhs.type, env))
        # Use LetRefNode to avoid side effects.
        let_ref_nodes.reverse()
        for t in let_ref_nodes:
            node = LetNode(t, node)
        return node

    def visit_ExprNode(self, node):
        # In-place assignments can't happen within an expression.
        return node


class AdjustDefByDirectives(CythonTransform, SkipDeclarations):
    """
    Adjust function and class definitions by the decorator directives:

    @cython.cfunc
    @cython.cclass
    @cython.ccall
    @cython.inline
    @cython.nogil
    @cython.critical_section
    """
    # list of directives that cause conversion to cclass
    converts_to_cclass = ('cclass', 'total_ordering', 'dataclasses.dataclass')

    def visit_ModuleNode(self, node):
        self.directives = node.directives
        self.in_py_class = False
        self.visitchildren(node)
        return node

    def visit_CompilerDirectivesNode(self, node):
        old_directives = self.directives
        self.directives = node.directives
        self.visitchildren(node)
        self.directives = old_directives
        return node

    def visit_DefNode(self, node):
        modifiers = []
        if 'inline' in self.directives:
            modifiers.append('inline')
        nogil = self.directives.get('nogil')
        with_gil = self.directives.get('with_gil')
        except_val = self.directives.get('exceptval')
        has_explicit_exc_clause = False if except_val is None else True
        return_type_node = self.directives.get('returns')
        if return_type_node is None and self.directives['annotation_typing']:
            return_type_node = node.return_type_annotation
            # for Python annotations, prefer safe exception handling by default
            if return_type_node is not None and except_val is None:
                except_val = (None, True)  # except *
        elif except_val is None:
            # backward compatible default: no exception check, unless there's also a "@returns" declaration
            except_val = (None, True if return_type_node else False)
        if self.directives.get('c_compile_guard') and 'cfunc' not in self.directives:
            error(node.pos, "c_compile_guard only allowed on C functions")
        if 'ccall' in self.directives:
            if 'cfunc' in self.directives:
                error(node.pos, "cfunc and ccall directives cannot be combined")
            if with_gil:
                error(node.pos, "ccall functions cannot be declared 'with_gil'")
            node = node.as_cfunction(
                overridable=True, modifiers=modifiers, nogil=nogil,
                returns=return_type_node, except_val=except_val, has_explicit_exc_clause=has_explicit_exc_clause)
            return self.visit(node)
        if 'cfunc' in self.directives:
            if self.in_py_class:
                error(node.pos, "cfunc directive is not allowed here")
            else:
                node = node.as_cfunction(
                    overridable=False, modifiers=modifiers, nogil=nogil, with_gil=with_gil,
                    returns=return_type_node, except_val=except_val, has_explicit_exc_clause=has_explicit_exc_clause)
                return self.visit(node)
        if 'inline' in modifiers:
            error(node.pos, "Python functions cannot be declared 'inline'")
        if nogil:
            # TODO: turn this into a "with gil" declaration.
            error(node.pos, "Python functions cannot be declared 'nogil'")
        if with_gil:
            error(node.pos, "Python functions cannot be declared 'with_gil'")
        self.visit_FuncDefNode(node)
        return node

    def visit_FuncDefNode(self, node):
        if "critical_section" in self.directives:
            value = self.directives["critical_section"]
            if value is not None:
                error(node.pos, "critical_section decorator does not take arguments")
            new_body = Nodes.CriticalSectionStatNode(
                node.pos,
                args=[ExprNodes.FirstArgumentForCriticalSectionNode(node.pos, func_node=node)],
                body=node.body
            )
            node.body = new_body
        self.visitchildren(node)
        return node

    def visit_LambdaNode(self, node):
        # No directives should modify lambdas or generator expressions (and also nothing in them).
        return node

    def visit_PyClassDefNode(self, node):
        if any(directive in self.directives for directive in self.converts_to_cclass):
            node = node.as_cclass()
            return self.visit(node)
        else:
            old_in_pyclass = self.in_py_class
            self.in_py_class = True
            self.visitchildren(node)
            self.in_py_class = old_in_pyclass
            return node

    def visit_CClassDefNode(self, node):
        old_in_pyclass = self.in_py_class
        self.in_py_class = False
        self.visitchildren(node)
        self.in_py_class = old_in_pyclass
        return node


class AlignFunctionDefinitions(CythonTransform):
    """
    This class takes the signatures from a .pxd file and applies them to
    the def methods in a .py file.
    """

    def visit_ModuleNode(self, node):
        self.scope = node.scope
        self.visitchildren(node)
        return node

    def visit_PyClassDefNode(self, node):
        pxd_def = self.scope.lookup(node.name)
        if pxd_def:
            if pxd_def.is_cclass:
                return self.visit_CClassDefNode(node.as_cclass(), pxd_def)
            elif not pxd_def.scope or not pxd_def.scope.is_builtin_scope:
                error(node.pos, "'%s' redeclared" % node.name)
                if pxd_def.pos:
                    error(pxd_def.pos, "previous declaration here")
                return None
        return node

    def visit_CClassDefNode(self, node, pxd_def=None):
        if pxd_def is None:
            pxd_def = self.scope.lookup(node.class_name)
        if pxd_def:
            if not pxd_def.defined_in_pxd:
                return node
            outer_scope = self.scope
            self.scope = pxd_def.type.scope
        self.visitchildren(node)
        if pxd_def:
            self.scope = outer_scope
        return node

    def visit_DefNode(self, node):
        pxd_def = self.scope.lookup(node.name)
        if pxd_def and (not pxd_def.scope or not pxd_def.scope.is_builtin_scope):
            if not pxd_def.is_cfunction:
                error(node.pos, "'%s' redeclared" % node.name)
                if pxd_def.pos:
                    error(pxd_def.pos, "previous declaration here")
                return None
            node = node.as_cfunction(pxd_def)
        # Enable this when nested cdef functions are allowed.
        # self.visitchildren(node)
        return node

    def visit_ExprNode(self, node):
        # ignore lambdas and everything else that appears in expressions
        return node


class AutoCpdefFunctionDefinitions(CythonTransform):

    def visit_ModuleNode(self, node):
        self.directives = node.directives
        self.imported_names = set()  # hack, see visit_FromImportStatNode()
        self.scope = node.scope
        self.visitchildren(node)
        return node

    def visit_DefNode(self, node):
        if (self.scope.is_module_scope and self.directives['auto_cpdef']
                and node.name not in self.imported_names
                and node.is_cdef_func_compatible()):
            # FIXME: cpdef-ing should be done in analyse_declarations()
            node = node.as_cfunction(scope=self.scope)
        return node

    def visit_CClassDefNode(self, node, pxd_def=None):
        if pxd_def is None:
            pxd_def = self.scope.lookup(node.class_name)
        if pxd_def:
            if not pxd_def.defined_in_pxd:
                return node
            outer_scope = self.scope
            self.scope = pxd_def.type.scope
        self.visitchildren(node)
        if pxd_def:
            self.scope = outer_scope
        return node

    def visit_FromImportStatNode(self, node):
        # hack to prevent conditional import fallback functions from
        # being cdpef-ed (global Python variables currently conflict
        # with imports)
        if self.scope.is_module_scope:
            for name, _ in node.items:
                self.imported_names.add(name)
        return node

    def visit_ExprNode(self, node):
        # ignore lambdas and everything else that appears in expressions
        return node


class RemoveUnreachableCode(CythonTransform):

    def visit_StatListNode(self, node):
        if not self.current_directives['remove_unreachable']:
            return node
        self.visitchildren(node)
        if len(node.stats) == 1 and isinstance(node.stats[0], Nodes.StatListNode) and not node.stats[0].stats:
            del node.stats[:]
        for idx, stat in enumerate(node.stats, 1):
            if stat.is_terminator:
                if idx < len(node.stats):
                    if self.current_directives['warn.unreachable']:
                        warning(node.stats[idx].pos, "Unreachable code", 2)
                    node.stats = node.stats[:idx]
                node.is_terminator = True
                break
        return node

    def visit_IfClauseNode(self, node):
        self.visitchildren(node)
        if node.body.is_terminator:
            node.is_terminator = True
        return node

    def visit_IfStatNode(self, node):
        self.visitchildren(node)
        if node.else_clause and node.else_clause.is_terminator:
            for clause in node.if_clauses:
                if not clause.is_terminator:
                    break
            else:
                node.is_terminator = True
        return node

    def visit_TryExceptStatNode(self, node):
        self.visitchildren(node)
        if node.body.is_terminator and node.else_clause:
            if self.current_directives['warn.unreachable']:
                warning(node.else_clause.pos, "Unreachable code", 2)
            node.else_clause = None
        return node

    def visit_TryFinallyStatNode(self, node):
        self.visitchildren(node)
        if node.finally_clause.is_terminator:
            node.is_terminator = True
        return node

    def visit_PassStatNode(self, node):
        """Eliminate useless PassStatNode"""
        # 'pass' statements often appear in a separate line and must be traced.
        if not self.current_directives['linetrace']:
            node = Nodes.StatListNode(pos=node.pos, stats=[])
        return node


class YieldNodeCollector(TreeVisitor):

    def __init__(self, excludes=[]):
        super().__init__()
        self.yields = []
        self.returns = []
        self.finallys = []
        self.excepts = []
        self.has_return_value = False
        self.has_yield = False
        self.has_await = False
        self.excludes = excludes

    def visit_Node(self, node):
        if node not in self.excludes:
            self.visitchildren(node)

    def visit_YieldExprNode(self, node):
        self.yields.append(node)
        self.has_yield = True
        self.visitchildren(node)

    def visit_AwaitExprNode(self, node):
        self.yields.append(node)
        self.has_await = True
        self.visitchildren(node)

    def visit_ReturnStatNode(self, node):
        self.visitchildren(node)
        if node.value:
            self.has_return_value = True
        self.returns.append(node)

    def visit_TryFinallyStatNode(self, node):
        self.visitchildren(node)
        self.finallys.append(node)

    def visit_TryExceptStatNode(self, node):
        self.visitchildren(node)
        self.excepts.append(node)

    def visit_ClassDefNode(self, node):
        pass

    def visit_FuncDefNode(self, node):
        pass

    def visit_LambdaNode(self, node):
        pass

    def visit_GeneratorExpressionNode(self, node):
        # node.loop iterator is evaluated outside the generator expression
        if isinstance(node.loop, Nodes._ForInStatNode):
            # Possibly should handle ForFromStatNode
            # but for now do nothing
            self.visit(node.loop.iterator)

    def visit_CArgDeclNode(self, node):
        # do not look into annotations
        # FIXME: support (yield) in default arguments (currently crashes)
        pass


class MarkClosureVisitor(CythonTransform):
    # In addition to marking closures this is also responsible to finding parts of the
    # generator iterable and marking them

    def visit_ModuleNode(self, node):
        self.needs_closure = False
        self.excludes = []
        self.visitchildren(node)
        return node

    def visit_FuncDefNode(self, node):
        self.needs_closure = False
        self.visitchildren(node)
        node.needs_closure = self.needs_closure
        self.needs_closure = True

        collector = YieldNodeCollector(self.excludes)
        collector.visitchildren(node)

        if node.is_async_def:
            coroutine_type = Nodes.AsyncDefNode
            if collector.has_yield:
                coroutine_type = Nodes.AsyncGenNode
                for yield_expr in collector.yields + collector.returns:
                    yield_expr.in_async_gen = True
            elif self.current_directives['iterable_coroutine']:
                coroutine_type = Nodes.IterableAsyncDefNode
        elif collector.has_await:
            found = next(y for y in collector.yields if y.is_await)
            error(found.pos, "'await' not allowed in generators (use 'yield')")
            return node
        elif collector.has_yield:
            coroutine_type = Nodes.GeneratorDefNode
        else:
            return node

        for i, yield_expr in enumerate(collector.yields, 1):
            yield_expr.label_num = i
        for retnode in collector.returns + collector.finallys + collector.excepts:
            retnode.in_generator = True

        gbody = Nodes.GeneratorBodyDefNode(
            pos=node.pos, name=node.name, body=node.body,
            is_coroutine_body=node.is_async_def,
            is_async_gen_body=node.is_async_def and collector.has_yield)
        coroutine = coroutine_type(
            pos=node.pos, name=node.name, args=node.args,
            star_arg=node.star_arg, starstar_arg=node.starstar_arg,
            doc=node.doc, decorators=node.decorators,
            gbody=gbody, lambda_name=node.lambda_name,
            return_type_annotation=node.return_type_annotation,
            is_generator_expression=node.is_generator_expression)
        return coroutine

    def visit_CFuncDefNode(self, node):
        self.needs_closure = False
        self.visitchildren(node)
        node.needs_closure = self.needs_closure
        self.needs_closure = True
        if node.needs_closure and node.overridable:
            error(node.pos, "closures inside cpdef functions not yet supported")
        return node

    def visit_LambdaNode(self, node):
        self.needs_closure = False
        self.visitchildren(node)
        node.needs_closure = self.needs_closure
        self.needs_closure = True
        return node

    def visit_ClassDefNode(self, node):
        self.visitchildren(node)
        self.needs_closure = True
        return node

    def visit_GeneratorExpressionNode(self, node):
        excludes = self.excludes
        if isinstance(node.loop, Nodes._ForInStatNode):
            self.excludes = [node.loop.iterator]
        node = self.visit_LambdaNode(node)
        self.excludes = excludes
        if not isinstance(node.loop, Nodes._ForInStatNode):
            # Possibly should handle ForFromStatNode
            # but for now do nothing
            return node
        itseq = node.loop.iterator.sequence
        # literals do not need replacing with an argument
        if itseq.is_literal:
            return node
        _GeneratorExpressionArgumentsMarker(node).visit(itseq)
        return node


class CreateClosureClasses(CythonTransform):
    # Output closure classes in module scope for all functions
    # that really need it.

    def __init__(self, context):
        super().__init__(context)
        self.path = []
        self.in_lambda = False

    def visit_ModuleNode(self, node):
        self.module_scope = node.scope
        self.visitchildren(node)
        return node

    def find_entries_used_in_closures(self, node):
        from_closure = []
        in_closure = []
        for scope in node.local_scope.iter_local_scopes():
            for name, entry in scope.entries.items():
                if not name:
                    continue
                if entry.from_closure:
                    from_closure.append((name, entry))
                elif entry.in_closure:
                    in_closure.append((name, entry))
        return from_closure, in_closure

    def create_class_from_scope(self, node, target_module_scope, inner_node=None):
        # move local variables into closure
        if node.is_generator:
            for scope in node.local_scope.iter_local_scopes():
                for entry in scope.entries.values():
                    if not (entry.from_closure or entry.is_pyglobal or entry.is_cglobal):
                        entry.in_closure = True

        from_closure, in_closure = self.find_entries_used_in_closures(node)
        in_closure.sort()

        # Now from the beginning
        node.needs_closure = False
        node.needs_outer_scope = False

        func_scope = node.local_scope
        cscope = node.entry.scope
        while cscope.is_py_class_scope or cscope.is_c_class_scope:
            cscope = cscope.outer_scope

        if not from_closure and (self.path or inner_node):
            if not inner_node:
                if not node.py_cfunc_node:
                    raise InternalError("DefNode does not have assignment node")
                inner_node = node.py_cfunc_node
            inner_node.needs_closure_code = False
            node.needs_outer_scope = False

        if node.is_generator:
            pass
        elif not in_closure and not from_closure:
            return
        elif not in_closure:
            func_scope.is_passthrough = True
            func_scope.scope_class = cscope.scope_class
            node.needs_outer_scope = True
            return

        # entry.cname can contain periods (eg. a derived C method of a class).
        # We want to use the cname as part of a C struct name, so we replace
        # periods with double underscores.
        as_name = '%s_%s' % (
            target_module_scope.next_id(Naming.closure_class_prefix),
            node.entry.cname.replace('.','__'))
        as_name = EncodedString(as_name)

        entry = target_module_scope.declare_c_class(
            name=as_name, pos=node.pos, defining=True,
            implementing=True)
        entry.type.is_final_type = True

        func_scope.scope_class = entry
        class_scope = entry.type.scope
        class_scope.is_internal = True
        class_scope.is_closure_class_scope = True
        if node.is_async_def or node.is_generator:
            # Generators need their closure intact during cleanup as they resume to handle GeneratorExit
            class_scope.directives['no_gc_clear'] = True
        if Options.closure_freelist_size:
            class_scope.directives['freelist'] = Options.closure_freelist_size

        if from_closure:
            assert cscope.is_closure_scope
            class_scope.declare_var(pos=node.pos,
                                    name=Naming.outer_scope_cname,
                                    cname=Naming.outer_scope_cname,
                                    type=cscope.scope_class.type,
                                    is_cdef=True)
            node.needs_outer_scope = True
        for name, entry in in_closure:
            closure_entry = class_scope.declare_var(
                pos=entry.pos,
                name=entry.name if not entry.in_subscope else None,
                cname=entry.cname,
                type=entry.type,
                is_cdef=True)
            if entry.is_declared_generic:
                closure_entry.is_declared_generic = 1
        node.needs_closure = True
        # Do it here because other classes are already checked
        target_module_scope.check_c_class(func_scope.scope_class)

    def visit_LambdaNode(self, node):
        if not isinstance(node.def_node, Nodes.DefNode):
            # fused function, an error has been previously issued
            return node

        was_in_lambda = self.in_lambda
        self.in_lambda = True
        self.create_class_from_scope(node.def_node, self.module_scope, node)
        self.visitchildren(node)
        self.in_lambda = was_in_lambda
        return node

    def visit_FuncDefNode(self, node):
        if self.in_lambda:
            self.visitchildren(node)
            return node
        if node.needs_closure or self.path:
            self.create_class_from_scope(node, self.module_scope)
            self.path.append(node)
            self.visitchildren(node)
            self.path.pop()
        return node

    def visit_GeneratorBodyDefNode(self, node):
        self.visitchildren(node)
        return node

    def visit_CFuncDefNode(self, node):
        if not node.overridable:
            return self.visit_FuncDefNode(node)
        else:
            self.visitchildren(node)
            return node

    def visit_GeneratorExpressionNode(self, node):
        node = _HandleGeneratorArguments()(node)
        return self.visit_LambdaNode(node)


class InjectGilHandling(VisitorTransform, SkipDeclarations):
    """
    Allow certain Python operations inside of nogil blocks by implicitly acquiring the GIL.

    Must run before the AnalyseDeclarationsTransform to make sure the GILStatNodes get
    set up, parallel sections know that the GIL is acquired inside of them, etc.
    """
    nogil = False

    # special node handling

    def _inject_gil_in_nogil(self, node):
        """Allow the (Python statement) node in nogil sections by wrapping it in a 'with gil' block."""
        if self.nogil:
            node = Nodes.GILStatNode(node.pos, state='gil', body=node)
        return node

    visit_RaiseStatNode = _inject_gil_in_nogil
    visit_PrintStatNode = _inject_gil_in_nogil  # sadly, not the function

    # further candidates:
    # def visit_ReraiseStatNode(self, node):

    # nogil tracking

    def visit_GILStatNode(self, node):
        was_nogil = self.nogil
        self.nogil = (node.state == 'nogil')
        self.visitchildren(node)
        self.nogil = was_nogil
        return node

    def visit_CFuncDefNode(self, node):
        was_nogil = self.nogil
        if isinstance(node.declarator, Nodes.CFuncDeclaratorNode):
            self.nogil = node.declarator.nogil and not node.declarator.with_gil
        self.visitchildren(node)
        self.nogil = was_nogil
        return node

    def visit_ParallelRangeNode(self, node):
        was_nogil = self.nogil
        self.nogil = node.nogil
        self.visitchildren(node)
        self.nogil = was_nogil
        return node

    def visit_ExprNode(self, node):
        # No special GIL handling inside of expressions for now.
        return node

    visit_Node = VisitorTransform.recurse_to_children


class GilCheck(VisitorTransform):
    """
    Call `node.gil_check(env)` on each node to make sure we hold the
    GIL when we need it.  Raise an error when on Python operations
    inside a `nogil` environment.

    Additionally, raise exceptions for closely nested with gil or with nogil
    statements. The latter would abort Python.
    """

    def __call__(self, root):
        self.env_stack = [root.scope]
        self.nogil_state = Nodes.NoGilState.HasGil

        self.nogil_state_at_current_gilstatnode = Nodes.NoGilState.HasGil
        return super().__call__(root)

    def _visit_scoped_children(self, node, nogil_state):
        was_nogil = self.nogil_state
        outer_attrs = node.outer_attrs
        if outer_attrs and len(self.env_stack) > 1:
            self.nogil_state = (
                Nodes.NoGilState.NoGil if self.env_stack[-2].nogil else Nodes.NoGilState.HasGil)
            self.visitchildren(node, outer_attrs)

        self.nogil_state = nogil_state
        self.visitchildren(node, attrs=None, exclude=outer_attrs)
        self.nogil_state = was_nogil

    def visit_FuncDefNode(self, node):
        self.env_stack.append(node.local_scope)
        inner_nogil = node.local_scope.nogil

        nogil_state = self.nogil_state
        if inner_nogil:
            self.nogil_state = Nodes.NoGilState.NoGilScope

        if inner_nogil and node.nogil_check:
            node.nogil_check(node.local_scope)

        self._visit_scoped_children(node, self.nogil_state)

        # FuncDefNodes can be nested, because a cpdef function contains a def function
        # inside it. Therefore restore to previous state
        self.nogil_state = nogil_state

        self.env_stack.pop()
        return node

    def visit_GILStatNode(self, node):
        if node.condition is not None:
            error(node.condition.pos,
                  "Non-constant condition in a "
                  "`with %s(<condition>)` statement" % node.state)
            return node

        if self.nogil_state and node.nogil_check:
            node.nogil_check()

        was_nogil = self.nogil_state
        is_nogil = (node.state == 'nogil')

        if was_nogil == is_nogil and not self.nogil_state == Nodes.NoGilState.NoGilScope:
            if not was_nogil:
                error(node.pos, "Trying to acquire the GIL while it is "
                                "already held.")
            else:
                error(node.pos, "Trying to release the GIL while it was "
                                "previously released.")
        if self.nogil_state == Nodes.NoGilState.NoGilScope:
            node.scope_gil_state_known = False

        if isinstance(node.finally_clause, Nodes.StatListNode):
            # The finally clause of the GILStatNode is a GILExitNode,
            # which is wrapped in a StatListNode. Just unpack that.
            node.finally_clause, = node.finally_clause.stats

        nogil_state_at_current_gilstatnode = self.nogil_state_at_current_gilstatnode
        self.nogil_state_at_current_gilstatnode = self.nogil_state
        nogil_state = Nodes.NoGilState.NoGil if is_nogil else Nodes.NoGilState.HasGil
        self._visit_scoped_children(node, nogil_state)
        self.nogil_state_at_current_gilstatnode = nogil_state_at_current_gilstatnode
        return node

    def visit_ParallelRangeNode(self, node):
        if node.nogil or self.nogil_state == Nodes.NoGilState.NoGilScope:
            node_was_nogil, node.nogil = node.nogil, False
            node = Nodes.GILStatNode(node.pos, state='nogil', body=node)
            if not node_was_nogil and self.nogil_state == Nodes.NoGilState.NoGilScope:
                # We're in a "nogil" function, but that doesn't prove we
                # didn't have the gil
                node.scope_gil_state_known = False
            return self.visit_GILStatNode(node)

        if not self.nogil_state:
            error(node.pos, "prange() can only be used without the GIL")
            # Forget about any GIL-related errors that may occur in the body
            return None

        node.nogil_check(self.env_stack[-1])
        self.visitchildren(node)
        return node

    def visit_ParallelWithBlockNode(self, node):
        if not self.nogil_state:
            error(node.pos, "The parallel section may only be used without "
                            "the GIL")
            return None
        if self.nogil_state == Nodes.NoGilState.NoGilScope:
            # We're in a "nogil" function but that doesn't prove we didn't
            # have the gil, so release it
            node = Nodes.GILStatNode(node.pos, state='nogil', body=node)
            node.scope_gil_state_known = False
            return self.visit_GILStatNode(node)

        if node.nogil_check:
            # It does not currently implement this, but test for it anyway to
            # avoid potential future surprises
            node.nogil_check(self.env_stack[-1])

        self.visitchildren(node)
        return node

    def visit_TryFinallyStatNode(self, node):
        """
        Take care of try/finally statements in nogil code sections.
        """
        if not self.nogil_state:
            return self.visit_Node(node)

        node.nogil_check = None
        node.is_try_finally_in_nogil = True
        self.visitchildren(node)
        return node

    def visit_CriticalSectionStatNode(self, node):
        # skip normal "try/finally node" handling
        return self.visit_Node(node)

    def visit_CythonLockStatNode(self, node):
        # skip normal "try/finally node" handling
        return self.visit_Node(node)

    def visit_GILExitNode(self, node):
        if self.nogil_state_at_current_gilstatnode == Nodes.NoGilState.NoGilScope:
            node.scope_gil_state_known = False
        self.visitchildren(node)
        return node

    def visit_Node(self, node):
        if self.env_stack and self.nogil_state and node.nogil_check:
            node.nogil_check(self.env_stack[-1])
        if node.outer_attrs:
            self._visit_scoped_children(node, self.nogil_state)
        else:
            self.visitchildren(node)
        if self.nogil_state:
            node.in_nogil_context = self.nogil_state
        return node

    def visit_SimpleCallNode(self, node):
        if (node.self and node.self.type.is_cython_lock_type and
                node.function.is_attribute and node.function.attribute == "acquire" and
                len(node.args) == 1):
            # For the cython lock types we can optimize if we know the GIL state.
            # (Remove this in the distant future when it's all PyMutexes because for these
            # it doesn't matter)
            suffix = None
            if self.nogil_state == Nodes.NoGilState.NoGil:
                suffix = "Nogil"
            elif self.nogil_state == Nodes.NoGilState.HasGil:
                suffix = "Gil"
            if suffix:
                node = ExprNodes.PythonCapiCallNode(
                    node.pos,
                    node.function.entry.cname + suffix,
                    node.function.type,
                    args=[node.self],
                )
        return self.visit_Node(node)


class CoerceCppTemps(EnvTransform, SkipDeclarations):
    """
    For temporary expression that are implemented using std::optional it's necessary the temps are
    assigned using `__pyx_t_x = value;` but accessed using `something = (*__pyx_t_x)`. This transform
    inserts a coercion node to take care of this, and runs absolutely last (once nothing else can be
    inserted into the tree)

    TODO: a possible alternative would be to split ExprNode.result() into ExprNode.rhs_result() and ExprNode.lhs_result()???
    """
    def visit_ModuleNode(self, node):
        if self.current_env().cpp:
            # skipping this makes it essentially free for C files
            self.visitchildren(node)
        return node

    def visit_ExprNode(self, node):
        self.visitchildren(node)
        if (self.current_env().directives['cpp_locals'] and
                node.result_in_temp() and node.type.is_cpp_class and
                # Fake references are not replaced with "std::optional()".
                not node.type.is_fake_reference):
            node = ExprNodes.CppOptionalTempCoercion(node)

        return node

    def visit_CppOptionalTempCoercion(self, node):
        return node

    def visit_CppIteratorNode(self, node):
        return node

    def visit_ExprStatNode(self, node):
        # Deliberately skip `expr` in ExprStatNode - we don't need to access it.
        self.visitchildren(node.expr)
        return node


class TransformBuiltinMethods(EnvTransform):
    """
    Replace Cython's own cython.* builtins by the corresponding tree nodes.
    Also handle some Python special builtin functions (e.g. super()/locals())
    that require introspection by the compiler.
    """
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.def_node_body_insertions = {}

    def visit_SingleAssignmentNode(self, node):
        if node.declaration_only:
            return None
        else:
            self.visitchildren(node)
            return node

    def visit_AttributeNode(self, node):
        self.visitchildren(node)
        return self.visit_cython_attribute(node)

    def visit_NameNode(self, node):
        if node.name == u'__class__':
            lenv = self.current_env()
            entry = lenv.lookup_here(u'__class__')
            if not entry:
                node = self._inject_class(node)
        return self.visit_cython_attribute(node)

    def visit_cython_attribute(self, node):
        attribute = node.as_cython_attribute()
        if attribute:
            if attribute == '__version__':
                from .. import __version__ as version
                node = ExprNodes.UnicodeNode(node.pos, value=EncodedString(version))
            elif attribute == 'NULL':
                node = ExprNodes.NullNode(node.pos)
            elif attribute in ('set', 'frozenset', 'staticmethod'):
                node = ExprNodes.NameNode(node.pos, name=EncodedString(attribute),
                                          entry=self.current_env().builtin_scope().lookup_here(attribute))
            elif PyrexTypes.parse_basic_type(attribute):
                pass
            elif self.context.cython_scope.lookup_qualified_name(attribute):
                pass
            else:
                error(node.pos, "'%s' not a valid cython attribute or is being used incorrectly" % attribute)
        return node

    def visit_ExecStatNode(self, node):
        lenv = self.current_env()
        self.visitchildren(node)
        if len(node.args) == 1:
            node.args.append(ExprNodes.GlobalsExprNode(node.pos))
            if not lenv.is_module_scope:
                node.args.append(
                    ExprNodes.LocalsExprNode(
                        node.pos, self.current_scope_node(), lenv))
        return node

    def _inject_locals(self, node, func_name):
        # locals()/dir()/vars() builtins
        lenv = self.current_env()
        entry = lenv.lookup_here(func_name)
        if entry:
            # not the builtin
            return node
        pos = node.pos
        if func_name in ('locals', 'vars'):
            if func_name == 'locals' and len(node.args) > 0:
                error(self.pos, "Builtin 'locals()' called with wrong number of args, expected 0, got %d"
                      % len(node.args))
                return node
            elif func_name == 'vars':
                if len(node.args) > 1:
                    error(self.pos, "Builtin 'vars()' called with wrong number of args, expected 0-1, got %d"
                          % len(node.args))
                if len(node.args) > 0:
                    return node  # nothing to do
            return ExprNodes.LocalsExprNode(pos, self.current_scope_node(), lenv)
        else:  # dir()
            if len(node.args) > 1:
                error(self.pos, "Builtin 'dir()' called with wrong number of args, expected 0-1, got %d"
                      % len(node.args))
            if len(node.args) > 0:
                # optimised in Builtin.py
                return node
            if lenv.is_py_class_scope or lenv.is_module_scope:
                if lenv.is_py_class_scope:
                    pyclass = self.current_scope_node()
                    locals_dict = ExprNodes.CloneNode(pyclass.dict)
                else:
                    locals_dict = ExprNodes.GlobalsExprNode(pos)
                return ExprNodes.SortedDictKeysNode(locals_dict)
            local_names = sorted(var.name for var in lenv.entries.values() if var.name)
            items = [ExprNodes.IdentifierStringNode(pos, value=var)
                     for var in local_names]
            return ExprNodes.ListNode(pos, args=items)

    def visit_PrimaryCmpNode(self, node):
        # special case: for in/not-in test, we do not need to sort locals()
        self.visitchildren(node)
        if node.operator in 'not_in':  # in/not_in
            if isinstance(node.operand2, ExprNodes.SortedDictKeysNode):
                arg = node.operand2.arg
                if isinstance(arg, ExprNodes.NoneCheckNode):
                    arg = arg.arg
                node.operand2 = arg
        return node

    def visit_CascadedCmpNode(self, node):
        return self.visit_PrimaryCmpNode(node)

    def _inject_eval(self, node, func_name):
        lenv = self.current_env()
        entry = lenv.lookup(func_name)
        if len(node.args) != 1 or (entry and not entry.is_builtin):
            return node
        # Inject globals and locals
        node.args.append(ExprNodes.GlobalsExprNode(node.pos))
        if not lenv.is_module_scope:
            node.args.append(
                ExprNodes.LocalsExprNode(
                    node.pos, self.current_scope_node(), lenv))
        return node

    def _inject_class(self, node):
        # bare __class__ reference inside function
        current_def_node = self.current_scope_node()

        if not isinstance(current_def_node, Nodes.FuncDefNode):
            return node

        # Go up the stack, find the first class node and its direct method (i.e. function).
        fdef_node = class_node = generator_node = None
        for stack_node, stack_scope in reversed(self.env_stack):
            if isinstance(stack_node, Nodes.ClassDefNode):
                class_node = stack_node
                class_scope = stack_scope
                break
            elif isinstance(stack_node, Nodes.GeneratorDefNode):
                generator_node = stack_node
                fdef_node = stack_node.gbody
                fdef_scope = stack_scope
            elif isinstance(stack_node, Nodes.FuncDefNode):
                fdef_node = stack_node
                fdef_scope = stack_scope

        if not fdef_node or not class_node:
            # failed to find a class or function
            return node

        # now we arrange to inject:
        #  __class__ = ... at the start of the def_node body
        # The advantage of doing it like this is that it automatically appears in locals()
        # and it can be captured by inner functions
        if fdef_node not in self.def_node_body_insertions:
            pos = fdef_node.body.pos
            if class_scope.is_c_class_scope:
                # c-classes can be resolved at compile-time, so they have a simpler
                # implementation
                rhs = ExprNodes.NameNode(
                            pos, name=class_node.scope.name,
                            entry=class_node.entry)
            elif class_scope.is_py_class_scope:
                rhs = ExprNodes.ClassCellNode(pos, is_generator=generator_node is not None)
                if generator_node:
                    generator_node.requires_classobj = True
                else:
                    fdef_node.requires_classobj = True
                class_node.class_cell.is_active = True
            else:
                return node  # should never happen

            assign_node = Nodes.SingleAssignmentNode(pos,
                lhs=ExprNodes.NameNode(pos, name=EncodedString("__class__")),
                rhs=rhs)

            assign_node.analyse_declarations(fdef_scope)

            assert fdef_node not in self.def_node_body_insertions
            self.def_node_body_insertions[fdef_node] = assign_node

        return node

    def _inject_super(self, node, func_name):
        lenv = self.current_env()
        entry = lenv.lookup_here(func_name)
        if entry or node.args:
            return node
        # Inject no-args super
        def_node = self.current_scope_node()
        if not isinstance(def_node, Nodes.DefNode) or not def_node.args or len(self.env_stack) < 2:
            return node
        class_node, class_scope = self.env_stack[-2]
        if class_scope.is_py_class_scope:
            def_node.requires_classobj = True
            class_node.class_cell.is_active = True
            node.args = [
                ExprNodes.ClassCellNode(
                    node.pos, is_generator=def_node.is_generator),
                ExprNodes.NameNode(node.pos, name=def_node.args[0].name)
                ]
        elif class_scope.is_c_class_scope:
            node.args = [
                ExprNodes.NameNode(
                    node.pos, name=class_node.scope.name,
                    entry=class_node.entry),
                ExprNodes.NameNode(node.pos, name=def_node.args[0].name)
                ]
        return node

    def _do_body_insertion(self, node):
        body_insertion = self.def_node_body_insertions.pop(node, None)
        if body_insertion:
            if isinstance(node.body, Nodes.StatListNode):
                node.body.stats.insert(0, body_insertion)
            else:
                node.body = Nodes.StatListNode(node.body.pos,
                                               stats=[body_insertion, node.body])

    def visit_FuncDefNode(self, node):
        node = super().visit_FuncDefNode(node)
        self._do_body_insertion(node)
        return node

    def visit_GeneratorBodyDefNode(self, node):
        node = super().visit_GeneratorBodyDefNode(node)
        self._do_body_insertion(node)
        return node

    def visit_SimpleCallNode(self, node):
        # cython.foo
        function = node.function.as_cython_attribute()
        if function:
            if function in InterpretCompilerDirectives.unop_method_nodes:
                if len(node.args) != 1:
                    error(node.function.pos, "%s() takes exactly one argument" % function)
                else:
                    node = InterpretCompilerDirectives.unop_method_nodes[function](
                        node.function.pos, operand=node.args[0])
            elif function in InterpretCompilerDirectives.binop_method_nodes:
                if len(node.args) != 2:
                    error(node.function.pos, "%s() takes exactly two arguments" % function)
                else:
                    node = InterpretCompilerDirectives.binop_method_nodes[function](
                        node.function.pos, operand1=node.args[0], operand2=node.args[1])
            elif function == 'cast':
                if len(node.args) != 2:
                    error(node.function.pos,
                          "cast() takes exactly two arguments and an optional typecheck keyword")
                else:
                    type = node.args[0].analyse_as_type(self.current_env())
                    if type:
                        node = ExprNodes.TypecastNode(
                            node.function.pos, type=type, operand=node.args[1], typecheck=False)
                    else:
                        error(node.args[0].pos, "Not a type")
            elif function == 'sizeof':
                if len(node.args) != 1:
                    error(node.function.pos, "sizeof() takes exactly one argument")
                else:
                    type = node.args[0].analyse_as_type(self.current_env())
                    if type:
                        node = ExprNodes.SizeofTypeNode(node.function.pos, arg_type=type)
                    else:
                        node = ExprNodes.SizeofVarNode(node.function.pos, operand=node.args[0])
            elif function == 'cmod':
                if len(node.args) != 2:
                    error(node.function.pos, "cmod() takes exactly two arguments")
                else:
                    node = ExprNodes.binop_node(node.function.pos, '%', node.args[0], node.args[1])
                    node.cdivision = True
            elif function == 'cdiv':
                if len(node.args) != 2:
                    error(node.function.pos, "cdiv() takes exactly two arguments")
                else:
                    node = ExprNodes.binop_node(node.function.pos, '/', node.args[0], node.args[1])
                    node.cdivision = True
            elif function == 'set':
                node.function = ExprNodes.NameNode(node.pos, name=EncodedString('set'))
            elif function == 'staticmethod':
                node.function = ExprNodes.NameNode(node.pos, name=EncodedString('staticmethod'))
            elif self.context.cython_scope.lookup_qualified_name(function):
                pass
            else:
                error(node.function.pos,
                      "'%s' not a valid cython language construct" % function)

        self.visitchildren(node)

        if isinstance(node, ExprNodes.SimpleCallNode) and node.function.is_name:
            func_name = node.function.name
            if func_name in ('dir', 'locals', 'vars'):
                return self._inject_locals(node, func_name)
            if func_name == 'eval':
                return self._inject_eval(node, func_name)
            if func_name == 'super':
                return self._inject_super(node, func_name)
        return node

    def visit_GeneralCallNode(self, node):
        function = node.function.as_cython_attribute()
        if function == 'cast':
            # NOTE: assuming simple tuple/dict nodes for positional_args and keyword_args
            args = node.positional_args.args
            kwargs = node.keyword_args.compile_time_value(None)
            if (len(args) != 2 or len(kwargs) > 1 or
                    (len(kwargs) == 1 and 'typecheck' not in kwargs)):
                error(node.function.pos,
                      "cast() takes exactly two arguments and an optional typecheck keyword")
            else:
                type = args[0].analyse_as_type(self.current_env())
                if type:
                    typecheck = kwargs.get('typecheck', False)
                    node = ExprNodes.TypecastNode(
                        node.function.pos, type=type, operand=args[1], typecheck=typecheck)
                else:
                    error(args[0].pos, "Not a type")

        self.visitchildren(node)
        return node


class ReplaceFusedTypeChecks(VisitorTransform):
    """
    This is not a transform in the pipeline. It is invoked on the specific
    versions of a cdef function with fused argument types. It filters out any
    type branches that don't match. e.g.

        if fused_t is mytype:
            ...
        elif fused_t in other_fused_type:
            ...
    """
    def __init__(self, local_scope):
        super().__init__()
        self.local_scope = local_scope
        # defer the import until now to avoid circular import time dependencies
        from .Optimize import ConstantFolding
        self.transform = ConstantFolding(reevaluate=True)

    def visit_IfStatNode(self, node):
        """
        Filters out any if clauses with false compile time type check
        expression.
        """
        self.visitchildren(node)
        return self.transform(node)

    def visit_GILStatNode(self, node):
        """
        Fold constant condition of GILStatNode.
        """
        self.visitchildren(node)
        return self.transform(node)

    def visit_PrimaryCmpNode(self, node):
        with Errors.local_errors(ignore=True):
            type1 = node.operand1.analyse_as_type(self.local_scope)
            type2 = node.operand2.analyse_as_type(self.local_scope)

        if type1 and type2:
            false_node = ExprNodes.BoolNode(node.pos, value=False)
            true_node = ExprNodes.BoolNode(node.pos, value=True)

            type1 = self.specialize_type(type1, node.operand1.pos)
            op = node.operator

            if op in ('is', 'is_not', '==', '!='):
                type2 = self.specialize_type(type2, node.operand2.pos)

                is_same = type1.same_as(type2)
                eq = op in ('is', '==')

                if (is_same and eq) or (not is_same and not eq):
                    return true_node

            elif op in ('in', 'not_in'):
                # We have to do an instance check directly, as operand2
                # needs to be a fused type and not a type with a subtype
                # that is fused. First unpack the typedef
                if isinstance(type2, PyrexTypes.CTypedefType):
                    type2 = type2.typedef_base_type

                if type1.is_fused:
                    error(node.operand1.pos, "Type is fused")
                elif not type2.is_fused:
                    error(node.operand2.pos,
                          "Can only use 'in' or 'not in' on a fused type")
                else:
                    types = PyrexTypes.get_specialized_types(type2)

                    for specialized_type in types:
                        if type1.same_as(specialized_type):
                            if op == 'in':
                                return true_node
                            else:
                                return false_node

                    if op == 'not_in':
                        return true_node

            return false_node

        return node

    def specialize_type(self, type, pos):
        try:
            return type.specialize(self.local_scope.fused_to_specific)
        except KeyError:
            error(pos, "Type is not specific")
            return type

    def visit_Node(self, node):
        self.visitchildren(node)
        return node


class DebugTransform(CythonTransform):
    """
    Write debug information for this Cython module.
    """

    def __init__(self, context, options, result):
        super().__init__(context)
        self.visited = set()
        # our treebuilder and debug output writer
        # (see Cython.Debugger.debug_output.CythonDebugWriter)
        self.tb = self.context.gdb_debug_outputwriter
        #self.c_output_file = options.output_file
        self.c_output_file = result.c_file

        # Closure support, basically treat nested functions as if the AST were
        # never nested
        self.nested_funcdefs = []

        # tells visit_NameNode whether it should register step-into functions
        self.register_stepinto = False

    def visit_ModuleNode(self, node):
        self.tb.module_name = node.full_module_name
        attrs = dict(
            module_name=node.full_module_name,
            filename=node.pos[0].filename,
            c_filename=self.c_output_file)

        self.tb.start('Module', attrs)

        # serialize functions
        self.tb.start('Functions')
        # First, serialize functions normally...
        self.visitchildren(node)

        # ... then, serialize nested functions
        for nested_funcdef in self.nested_funcdefs:
            self.visit_FuncDefNode(nested_funcdef)

        self.register_stepinto = True
        self.serialize_modulenode_as_function(node)
        self.register_stepinto = False
        self.tb.end('Functions')

        # 2.3 compatibility. Serialize global variables
        self.tb.start('Globals')
        entries = {}

        for k, v in node.scope.entries.items():
            if (v.qualified_name not in self.visited and not
                    v.name.startswith('__pyx_') and not
                    v.type.is_cfunction and not
                    v.type.is_extension_type):
                entries[k]= v

        self.serialize_local_variables(entries)
        self.tb.end('Globals')
        # self.tb.end('Module') # end Module after the line number mapping in
        # Cython.Compiler.ModuleNode.ModuleNode._serialize_lineno_map
        return node

    def visit_FuncDefNode(self, node):
        self.visited.add(node.local_scope.qualified_name)

        if getattr(node, 'is_wrapper', False):
            return node

        if self.register_stepinto:
            self.nested_funcdefs.append(node)
            return node

        # node.entry.visibility = 'extern'
        if node.py_func is None:
            pf_cname = ''
        else:
            pf_cname = node.py_func.entry.func_cname

        # For functions defined using def, cname will be pyfunc_cname=__pyx_pf_*
        # For functions defined using cpdef or cdef, cname will be func_cname=__pyx_f_*
        # In all cases, cname will be the name of the function containing the actual code
        cname = node.entry.pyfunc_cname or node.entry.func_cname

        attrs = dict(
            name=node.entry.name or getattr(node, 'name', '<unknown>'),
            cname=cname,
            pf_cname=pf_cname,
            qualified_name=node.local_scope.qualified_name,
            lineno=str(node.pos[1]))

        self.tb.start('Function', attrs=attrs)

        self.tb.start('Locals')
        self.serialize_local_variables(node.local_scope.entries)
        self.tb.end('Locals')

        self.tb.start('Arguments')
        for arg in node.local_scope.arg_entries:
            self.tb.start(arg.name)
            self.tb.end(arg.name)
        self.tb.end('Arguments')

        self.tb.start('StepIntoFunctions')
        self.register_stepinto = True
        self.visitchildren(node)
        self.register_stepinto = False
        self.tb.end('StepIntoFunctions')
        self.tb.end('Function')

        return node

    def visit_NameNode(self, node):
        if (self.register_stepinto and
                node.type is not None and
                node.type.is_cfunction and
                getattr(node, 'is_called', False) and
                node.entry.func_cname is not None):
            # don't check node.entry.in_cinclude, as 'cdef extern: ...'
            # declared functions are not 'in_cinclude'.
            # This means we will list called 'cdef' functions as
            # "step into functions", but this is not an issue as they will be
            # recognized as Cython functions anyway.
            attrs = dict(name=node.entry.func_cname)
            self.tb.start('StepIntoFunction', attrs=attrs)
            self.tb.end('StepIntoFunction')

        self.visitchildren(node)
        return node

    def serialize_modulenode_as_function(self, node):
        """
        Serialize the module-level code as a function so the debugger will know
        it's a "relevant frame" and it will know where to set the breakpoint
        for 'break modulename'.
        """
        self._serialize_modulenode_as_function(node, dict(
            name=node.full_module_name.rpartition('.')[-1],
            cname=node.module_init_func_cname(),
            pf_cname='',
            # Ignore the qualified_name, breakpoints should be set using
            # `cy break modulename:lineno` for module-level breakpoints.
            qualified_name='',
            lineno='1',
            is_initmodule_function="True",
        ))

    def _serialize_modulenode_as_function(self, node, attrs):
        self.tb.start('Function', attrs=attrs)

        self.tb.start('Locals')
        self.serialize_local_variables(node.scope.entries)
        self.tb.end('Locals')

        self.tb.start('Arguments')
        self.tb.end('Arguments')

        self.tb.start('StepIntoFunctions')
        self.register_stepinto = True
        self.visitchildren(node)
        self.register_stepinto = False
        self.tb.end('StepIntoFunctions')

        self.tb.end('Function')

    def serialize_local_variables(self, entries):
        for entry in entries.values():
            if not entry.cname:
                # not a local variable
                continue
            if entry.type.is_pyobject:
                vartype = 'PythonObject'
            else:
                vartype = 'CObject'

            if entry.from_closure:
                # We're dealing with a closure where a variable from an outer
                # scope is accessed, get it from the scope object.
                cname = '%s->%s' % (Naming.cur_scope_cname,
                                    entry.outer_entry.cname)

                qname = '%s.%s.%s' % (entry.scope.outer_scope.qualified_name,
                                      entry.scope.name,
                                      entry.name)
            elif entry.in_closure:
                cname = '%s->%s' % (Naming.cur_scope_cname,
                                    entry.cname)
                qname = entry.qualified_name
            else:
                cname = entry.cname
                qname = entry.qualified_name

            if not entry.pos:
                # this happens for variables that are not in the user's code,
                # e.g. for the global __builtins__, __doc__, etc. We can just
                # set the lineno to 0 for those.
                lineno = '0'
            else:
                lineno = str(entry.pos[1])

            attrs = dict(
                name=entry.name,
                cname=cname,
                qualified_name=qname,
                type=vartype,
                lineno=lineno)

            self.tb.start('LocalVar', attrs)
            self.tb.end('LocalVar')


class HasNoExceptionHandlingVisitor(TreeVisitor):
    """
    Used by finalExceptClauseNode to work out if the body
    needs to handle exceptions at all. This includes:

    1. Can raise an exception.
    2. May try to access the traceback.
    """
    def __init__(self):
        self.uses_no_exceptions = True
        self.assignment_lhs = None
        super().__init__()

    def __call__(self, node) -> bool:
        self.visit(node)
        return self.uses_no_exceptions

    def visit_Node(self, node):
        self.uses_no_exceptions = False  # In general, nodes use exceptions

    def visit_ExprStatNode(self, node):
        self.visitchildren(node)

    def visit_StatListNode(self, node):
        self.visitchildren(node)

    def visit_ExprNode(self, node):
        if not node.is_literal:
            self.uses_no_exceptions = False

    def visit_CallNode(self, node):
        # Implement this to make the behaviour as explicit as possible.
        # Even noexcept functions might end up printing a traceback.
        self.uses_no_exceptions = False

    def visit_PassStatNode(self, node):
        pass  # Does nothing.  Good.

    def visit_ReturnStatNode(self, node):
        if not self.uses_no_exceptions:
            return  # shortcut
        self.visitchildren(node)

    def visit_SingleAssignmentNode(self, node):
        if not self.uses_no_exceptions:
            return  # shortcut
        self.assignment_lhs = node.lhs
        self.visit(node.lhs)
        self.assignment_lhs = None
        rhs_type = node.rhs.type
        if not (rhs_type.is_numeric or rhs_type.is_pyobject or rhs_type.is_memoryviewslice):
            # Treat everything we haven't explicitly thought about as potentially dubious.
            # cpp classes may have non-trivial assignment operators for example.
            self.uses_no_exceptions = False
        if not self.uses_no_exceptions:
            return
        self.visitchildren(node, exclude=["lhs"])

    def visit_NameNode(self, node):
        if not self.uses_no_exceptions:
            return  # shortcut
        entry = node.entry
        if self.assignment_lhs is node:
            if not (entry.is_cglobal or entry.is_arg or
                    entry.is_local or entry.in_closure or entry.from_closure):
                self.uses_no_exceptions = False
                return
        else:
            if entry.is_cglobal:
                if entry.is_cpp_optional and node.initialized_check:
                    # Otherwise, reading C globals should be safe.
                    self.uses_no_exceptions = False
                    return
            elif entry.is_arg or entry.is_local or entry.in_closure or entry.from_closure:
                if (node.cf_is_null or node.cf_maybe_null) and not node.type.is_numeric:
                    # The logic here is slightly simpler than for NameNode error checking.
                    # This gives a few false negatives (which is always the safe thing to do)
                    # for memoryviews and cpp_optionals
                    self.uses_no_exceptions = False
                    return
            else:
                # Probably a py_global.
                self.uses_no_exceptions = False
                return

    def visit_AttributeNode(self, node):
        if node.is_py_attr:
            self.uses_no_exceptions = False
        elif (node.type.is_memoryviewslice or node.entry.is_cpp_optional) and self.assignment_lhs is not node:
            # Memoryviewslices and cpp_optional are OK as a target, but reading them involves checks.
            # (Although cpp optionals are currently banned elsewhere
            # because C++ classes may have non-trivial assignment).
            self.uses_no_exceptions = False
        # Python objects just need an incref and simple C types are fine, too. Others may not be.
        if not (node.type.is_pyobject or node.type.is_numeric or node.type.is_memoryviewslice):
            self.uses_no_exceptions = False
        if self.uses_no_exceptions:
            self.visitchildren(node)

    def visit_IndexNode(self, node):
        if not (node.base.type.is_array or node.base.type.is_ptr):
            self.uses_no_exceptions = False
        if not self.uses_no_exceptions:
            return
        self.visitchildren(node)

    def visit_CoerceToTempNode(self, node):
        self.visitchildren(node)
