""" Computes the Control Flow Graph of a function. """

from pythran.passmanager import FunctionAnalysis
from pythran.utils import isnum
from pythran.graph import DiGraph

import gast as ast


def is_true_predicate(node):
    # FIXME: there may be more patterns here
    if isnum(node) and node.value:
        return True
    if isinstance(node, ast.Attribute) and node.attr == 'True':
        return True
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)) and node.elts:
        return True
    if isinstance(node, ast.Dict) and node.keys:
        return True
    return False


class CFG(FunctionAnalysis):
    """
    Computes the Control Flow Graph of a function.

    The processing of a node yields a pair containing
    * the OUT nodes, to be linked with the IN nodes of the successor
    * the RAISE nodes, nodes that stop the control flow (exception/break/...)
    """

    #: The sink node in the control flow graph.
    #:
    #: The predecessors of this node are those AST nodes that terminate
    #: control flow without a return statement.
    NIL = object()

    ResultType = DiGraph

    def visit_FunctionDef(self, node):
        """OUT = node, RAISES = ()"""
        # the function itself is the entry point
        self.result.add_node(node)
        currs = (node,)
        for n in node.body:
            self.result.add_node(n)
            for curr in currs:
                self.result.add_edge(curr, n)
            currs, _ = self.visit(n)
        # add an edge to NIL for nodes that end the control flow
        # without a return
        self.result.add_node(CFG.NIL)
        for curr in currs:
            self.result.add_edge(curr, CFG.NIL)
        return (node,), ()

    def visit_Pass(self, node):
        """OUT = node, RAISES = ()"""
        return (node,), ()

    # All these nodes have the same behavior as pass
    visit_Assign = visit_AnnAssign = visit_AugAssign = visit_Import = visit_Pass
    visit_Expr = visit_Print = visit_ImportFrom = visit_Pass
    visit_Yield = visit_Delete = visit_Pass
    visit_Nonlocal = visit_Pass

    def visit_Return(self, node):
        """OUT = (), RAISES = ()"""
        return (), ()

    def visit_For(self, node):
        """
        OUT = (node,) + last body statements
        RAISES = body's that are not break or continue
        """
        currs = (node,)
        break_currs = tuple()
        raises = ()
        # handle body
        for n in node.body:
            self.result.add_node(n)
            for curr in currs:
                self.result.add_edge(curr, n)
            currs, nraises = self.visit(n)
            for nraise in nraises:
                if isinstance(nraise, ast.Break):
                    break_currs += (nraise,)
                elif isinstance(nraise, ast.Continue):
                    self.result.add_edge(nraise, node)
                else:
                    raises += (nraise,)
        # add the backward loop
        for curr in currs:
            self.result.add_edge(curr, node)

        # the else statement if needed
        if node.orelse:
            for n in node.orelse:
                self.result.add_node(n)
                for curr in currs:
                    self.result.add_edge(curr, n)
                currs, nraises = self.visit(n)
        else:
            currs = node,

        # while only
        if isinstance(node, ast.While):
            if is_true_predicate(node.test):
                return break_currs, raises
            else:
                return break_currs + currs, raises

        # for only
        return break_currs + currs, raises

    visit_While = visit_For

    def visit_If(self, node):
        """
        OUT = true branch U false branch
        RAISES = true branch U false branch
        """
        currs = (node,)
        raises = ()

        # true branch
        for n in node.body:
            self.result.add_node(n)
            for curr in currs:
                self.result.add_edge(curr, n)
            currs, nraises = self.visit(n)
            raises += nraises

        # false branch
        tcurrs = currs
        traises = raises
        currs = (node,)
        for n in node.orelse:
            self.result.add_node(n)
            for curr in currs:
                self.result.add_edge(curr, n)
            currs, nraises = self.visit(n)
            raises = traises + nraises

        if is_true_predicate(node.test):
            return tcurrs, raises

        return tcurrs + currs, raises

    def visit_Raise(self, node):
        """OUT = (), RAISES = (node)"""
        return (), (node,)

    visit_Break = visit_Continue = visit_Raise

    def visit_Assert(self, node):
        """OUT = RAISES = (node)"""
        return (node,), (node,)

    def visit_Try(self, node):
        """
        OUT = body's U handler's
        RAISES = handler's
        this equation is not has good has it could be...
        but we need type information to be more accurate
        """
        currs = (node,)
        raises = ()
        for handler in node.handlers:
            self.result.add_node(handler)
        for n in node.body:
            self.result.add_node(n)
            for curr in currs:
                self.result.add_edge(curr, n)
            currs, nraises = self.visit(n)
            for nraise in nraises:
                if isinstance(nraise, ast.Raise):
                    for handler in node.handlers:
                        self.result.add_edge(nraise, handler)
                else:
                    raises += (nraise,)
        for handler in node.handlers:
            ncurrs, nraises = self.visit(handler)
            currs += ncurrs
            raises += nraises
        return currs, raises

    def visit_ExceptHandler(self, node):
        """OUT = body's, RAISES = body's"""
        currs = (node,)
        raises = ()
        for n in node.body:
            self.result.add_node(n)
            for curr in currs:
                self.result.add_edge(curr, n)
            currs, nraises = self.visit(n)
            raises += nraises
        return currs, raises
