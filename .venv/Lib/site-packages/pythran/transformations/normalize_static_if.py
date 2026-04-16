""" NormalizeStaticIf adds support for static guards. """

from pythran.analyses import (ImportedIds, HasReturn, IsAssigned, CFG,
                              HasBreak, HasContinue, DefUseChains, Ancestors,
                              StaticExpressions, HasStaticExpression)
from pythran.passmanager import Transformation
from pythran.syntax import PythranSyntaxError

import gast as ast
from copy import deepcopy

LOOP_NONE, EARLY_RET, LOOP_BREAK, LOOP_CONT = range(4)


def outline(name, formal_parameters, out_parameters, stmts,
            has_return, has_break, has_cont):

    args = ast.arguments(
        [ast.Name(fp, ast.Param(), None, None) for fp in formal_parameters],
        [], None, [], [], None, [])

    if isinstance(stmts, ast.expr):
        assert not out_parameters, "no out parameters with expr"
        fdef = ast.FunctionDef(name, args, [ast.Return(stmts)], [], None, None)
    else:
        fdef = ast.FunctionDef(name, args, stmts, [], None, None)

        # this is part of a huge trick that plays with delayed type inference
        # it basically computes the return type based on out parameters, and
        # the return statement is unconditionally added so if we have other
        # returns, there will be a computation of the output type based on the
        # __combined of the regular return types and this one The original
        # returns have been patched above to have a different type that
        # cunningly combines with this output tuple
        #
        # This is the only trick I found to let pythran compute both the output
        # variable type and the early return type. But hey, a dirty one :-/

        stmts.append(
            ast.Return(
                ast.Tuple(
                    [ast.Name(fp, ast.Load(), None, None)
                     for fp in out_parameters],
                    ast.Load()
                )
            )
        )
        if has_return:
            pr = PatchReturn(stmts[-1], has_break or has_cont)
            pr.generic_visit(fdef)

        if has_break or has_cont:
            if not has_return:
                stmts[-1].value = ast.Tuple([ast.Constant(LOOP_NONE, None),
                                             stmts[-1].value],
                                            ast.Load())
            pbc = PatchBreakContinue(stmts[-1])
            pbc.generic_visit(fdef)

    return fdef


class PatchReturn(ast.NodeTransformer):

    def __init__(self, guard, has_break_or_cont):
        self.guard = guard
        self.has_break_or_cont = has_break_or_cont

    def visit_FunctionDef(self, node):
        return node

    def visit_Return(self, node):
        if node is self.guard:
            holder = "StaticIfNoReturn"
        else:
            holder = "StaticIfReturn"

        value = node.value

        return ast.Return(
            ast.Call(
                ast.Attribute(
                    ast.Attribute(
                        ast.Name("builtins", ast.Load(), None, None),
                        "pythran",
                        ast.Load()),
                    holder,
                    ast.Load()),
                [value] if value else [ast.Constant(None, None)],
                []))


class PatchBreakContinue(ast.NodeTransformer):

    def __init__(self, guard):
        self.guard = guard

    def visit_FunctionDef(self, node):
        return node

    def visit_For(self, node):
        return node

    def visit_While(self, node):
        return node

    def patch_Control(self, node, flag):
        new_node = deepcopy(self.guard)
        ret_val = new_node.value
        if isinstance(ret_val, ast.Call):
            if flag == LOOP_BREAK:
                ret_val.func.attr = "StaticIfBreak"
            else:
                ret_val.func.attr = "StaticIfCont"
        else:
            new_node.value.elts[0].value = flag
        return new_node

    def visit_Break(self, node):
        return self.patch_Control(node, LOOP_BREAK)

    def visit_Continue(self, node):
        return self.patch_Control(node, LOOP_CONT)


class NormalizeStaticIf(Transformation[StaticExpressions, Ancestors, DefUseChains]):

    def visit_Module(self, node):
        self.new_functions = []
        self.funcs = []
        self.cfgs = []
        self.generic_visit(node)
        node.body.extend(self.new_functions)
        return node

    def escaping_ids(self, scope_stmt, stmts):
        'gather sets of identifiers defined in stmts and used out of it'
        assigned_nodes = self.gather(IsAssigned, self.make_fake(stmts))
        escaping = set()
        for assigned_node in assigned_nodes:
            head = self.def_use_chains.chains[assigned_node]
            for user in head.users():
                if scope_stmt not in self.ancestors[user.node]:
                    escaping.add(head.name())
        return escaping

    @staticmethod
    def make_fake(stmts):
        return ast.If(ast.Constant(0, None), stmts, [])

    @staticmethod
    def make_dispatcher(static_expr, func_true, func_false,
                        imported_ids):
        dispatcher_args = [static_expr,
                           ast.Name(func_true.name, ast.Load(), None, None),
                           ast.Name(func_false.name, ast.Load(), None, None)]

        dispatcher = ast.Call(
            ast.Attribute(
                ast.Attribute(
                    ast.Name("builtins", ast.Load(), None, None),
                    "pythran",
                    ast.Load()),
                "static_if",
                ast.Load()),
            dispatcher_args, [])

        actual_call = ast.Call(
            dispatcher,
            [ast.Name(ii, ast.Load(), None, None) for ii in imported_ids],
            [])

        return actual_call

    def true_name(self):
        return "$isstatic{}".format(len(self.new_functions) + 0)

    def false_name(self):
        return "$isstatic{}".format(len(self.new_functions) + 1)

    def visit_FunctionDef(self, node):
        self.cfgs.append(self.gather(CFG, node))
        self.funcs.append(node)
        onode = self.generic_visit(node)
        self.funcs.pop()
        self.cfgs.pop()
        return onode

    def visit_IfExp(self, node):
        self.generic_visit(node)

        if node.test not in self.static_expressions:
            return node

        imported_ids = sorted(self.gather(ImportedIds, node))

        func_true = outline(self.true_name(), imported_ids, [],
                            node.body, False, False, False)
        func_false = outline(self.false_name(), imported_ids, [],
                             node.orelse, False, False, False)
        self.new_functions.extend((func_true, func_false))

        actual_call = self.make_dispatcher(node.test, func_true,
                                           func_false, imported_ids)

        return actual_call

    def make_control_flow_handlers(self, cont_n, status_n, expected_return,
                                   has_cont, has_break):
        '''
        Create the statements in charge of gathering control flow information
        for the static_if result, and executes the expected control flow
        instruction
        '''
        if expected_return:
            assign = cont_ass = [ast.Assign(
                [ast.Tuple(expected_return, ast.Store())],
                ast.Name(cont_n, ast.Load(), None, None), None)]
        else:
            assign = cont_ass = []

        if has_cont:
            cmpr = ast.Compare(ast.Name(status_n, ast.Load(), None, None),
                               [ast.Eq()], [ast.Constant(LOOP_CONT, None)])
            cont_ass = [ast.If(cmpr,
                               deepcopy(assign) + [ast.Continue()],
                               cont_ass)]
        if has_break:
            cmpr = ast.Compare(ast.Name(status_n, ast.Load(), None, None),
                               [ast.Eq()], [ast.Constant(LOOP_BREAK, None)])
            cont_ass = [ast.If(cmpr,
                               deepcopy(assign) + [ast.Break()],
                               cont_ass)]
        return cont_ass

    def visit_If(self, node):
        if node.test not in self.static_expressions:
            return self.generic_visit(node)

        imported_ids = self.gather(ImportedIds, node)

        assigned_ids_left = self.escaping_ids(node, node.body)
        assigned_ids_right = self.escaping_ids(node, node.orelse)
        assigned_ids_both = assigned_ids_left.union(assigned_ids_right)

        imported_ids.update(i for i in assigned_ids_left
                            if i not in assigned_ids_right)
        imported_ids.update(i for i in assigned_ids_right
                            if i not in assigned_ids_left)
        imported_ids = sorted(imported_ids)

        assigned_ids = sorted(assigned_ids_both)

        fbody = self.make_fake(node.body)
        true_has_return = self.gather(HasReturn, fbody)
        true_has_break = self.gather(HasBreak, fbody)
        true_has_cont = self.gather(HasContinue, fbody)

        felse = self.make_fake(node.orelse)
        false_has_return = self.gather(HasReturn, felse)
        false_has_break = self.gather(HasBreak, felse)
        false_has_cont = self.gather(HasContinue, felse)

        has_return = true_has_return or false_has_return
        has_break = true_has_break or false_has_break
        has_cont = true_has_cont or false_has_cont

        self.generic_visit(node)

        func_true = outline(self.true_name(), imported_ids, assigned_ids,
                            node.body, has_return, has_break, has_cont)
        func_false = outline(self.false_name(), imported_ids, assigned_ids,
                             node.orelse, has_return, has_break, has_cont)
        self.new_functions.extend((func_true, func_false))

        actual_call = self.make_dispatcher(node.test,
                                           func_true, func_false, imported_ids)

        # variable modified within the static_if
        expected_return = [ast.Name(ii, ast.Store(), None, None)
                           for ii in assigned_ids]

        self.update = True

        # name for various variables resulting from the static_if
        n = len(self.new_functions)
        status_n = "$status{}".format(n)
        return_n = "$return{}".format(n)
        cont_n = "$cont{}".format(n)

        if has_return:
            cfg = self.cfgs[-1]
            always_return = all(isinstance(x, (ast.Return, ast.Yield))
                                for x in cfg[node])
            always_return &= true_has_return and false_has_return

            fast_return = [ast.Name(status_n, ast.Store(), None, None),
                           ast.Name(return_n, ast.Store(), None, None),
                           ast.Name(cont_n, ast.Store(), None, None)]

            if always_return:
                return [ast.Assign([ast.Tuple(fast_return, ast.Store())],
                                   actual_call, None),
                        ast.Return(ast.Name(return_n, ast.Load(), None, None))]
            else:
                cont_ass = self.make_control_flow_handlers(cont_n, status_n,
                                                           expected_return,
                                                           has_cont, has_break)

                cmpr = ast.Compare(ast.Name(status_n, ast.Load(), None, None),
                                   [ast.Eq()], [ast.Constant(EARLY_RET, None)])
                return [ast.Assign([ast.Tuple(fast_return, ast.Store())],
                                   actual_call, None),
                        ast.If(cmpr,
                               [ast.Return(ast.Name(return_n, ast.Load(),
                                                    None, None))],
                               cont_ass)]
        elif has_break or has_cont:
            cont_ass = self.make_control_flow_handlers(cont_n, status_n,
                                                       expected_return,
                                                       has_cont, has_break)

            fast_return = [ast.Name(status_n, ast.Store(), None, None),
                           ast.Name(cont_n, ast.Store(), None, None)]
            return [ast.Assign([ast.Tuple(fast_return, ast.Store())],
                               actual_call, None)] + cont_ass
        elif expected_return:
            return ast.Assign([ast.Tuple(expected_return, ast.Store())],
                              actual_call, None)
        else:
            return ast.Expr(actual_call)


class SplitStaticExpression(Transformation[StaticExpressions]):

    def visit_Cond(self, node):
        '''
        generic expression splitting algorithm. Should work for ifexp and if
        using W(rap) and U(n)W(rap) to manage difference between expr and stmt

        The idea is to split a BinOp in three expressions:
            1. a (possibly empty) non-static expr
            2. an expr containing a static expr
            3. a (possibly empty) non-static expr
        Once split, the if body is refactored to keep the semantic,
        and then recursively split again, until all static expr are alone in a
        test condition
        '''
        NodeTy = type(node)
        if NodeTy is ast.IfExp:
            def W(x):
                return x

            def UW(x):
                return x
        else:
            def W(x):
                return [x]

            def UW(x):
                return x[0]

        has_static_expr = self.gather(HasStaticExpression, node.test)

        if not has_static_expr:
            return self.generic_visit(node)

        if node.test in self.static_expressions:
            return self.generic_visit(node)

        if not isinstance(node.test, ast.BinOp):
            return self.generic_visit(node)

        before, static = [], []
        values = [node.test.right, node.test.left]

        def has_static_expression(n):
            return self.gather(HasStaticExpression, n)

        while values and not has_static_expression(values[-1]):
            before.append(values.pop())

        while values and has_static_expression(values[-1]):
            static.append(values.pop())

        after = list(reversed(values))

        test_before = NodeTy(None, None, None)
        if before:
            assert len(before) == 1
            test_before.test = before[0]

        test_static = NodeTy(None, None, None)
        if static:
            test_static.test = static[0]
            if len(static) > 1:
                if after:
                    assert len(after) == 1
                    after = [ast.BinOp(static[1], node.test.op, after[0])]
                else:
                    after = static[1:]

        test_after = NodeTy(None, None, None)
        if after:
            assert len(after) == 1
            test_after.test = after[0]

        if isinstance(node.test.op, ast.BitAnd):
            if after:
                test_after.body = deepcopy(node.body)
                test_after.orelse = deepcopy(node.orelse)
                test_after = W(test_after)
            else:
                test_after = deepcopy(node.body)

            if static:
                test_static.body = test_after
                test_static.orelse = deepcopy(node.orelse)
                test_static = W(test_static)
            else:
                test_static = test_after

            if before:
                test_before.body = test_static
                test_before.orelse = node.orelse
                node = test_before
            else:
                node = UW(test_static)

        elif isinstance(node.test.op, ast.BitOr):
            if after:
                test_after.body = deepcopy(node.body)
                test_after.orelse = deepcopy(node.orelse)
                test_after = W(test_after)
            else:
                test_after = deepcopy(node.orelse)

            if static:
                test_static.body = deepcopy(node.body)
                test_static.orelse = test_after
                test_static = W(test_static)
            else:
                test_static = test_after

            if before:
                test_before.body = deepcopy(node.body)
                test_before.orelse = test_static
                node = test_before
            else:
                node = UW(test_static)
        else:
            raise PythranSyntaxError("operator not supported in a static if",
                                     node)

        self.update = True
        return self.generic_visit(node)

    visit_If = visit_IfExp = visit_Cond
