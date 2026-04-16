""" LazynessAnalysis returns number of time a name is use. """

from pythran.analyses.aliases import Aliases
from pythran.analyses.argument_effects import ArgumentEffects
from pythran.analyses.identifiers import Identifiers
from pythran.analyses.pure_expressions import PureExpressions
from pythran.passmanager import FunctionAnalysis
from pythran.syntax import PythranSyntaxError
from pythran.utils import get_variable, isattr
import pythran.metadata as md
import pythran.openmp as openmp

import gast as ast
import sys


class LazynessAnalysis(FunctionAnalysis[ArgumentEffects, Aliases, PureExpressions]):

    """
    Returns number of time a name is used.

    +inf if it is use in a
    loop, if a variable used to compute it is modify before
    its last use or if it is use in a function call (as it is not an
    interprocedural analysis)

    >>> import gast as ast, sys
    >>> from pythran import passmanager, backend
    >>> code = "def foo(): c = 1; a = c + 2; c = 2; b = c + c + a; return b"
    >>> node = ast.parse(code)
    >>> pm = passmanager.PassManager("test")
    >>> res = pm.gather(LazynessAnalysis, node)
    >>> res['a'], res['b'], res['c']
    (inf, 1, 2)
    >>> code = '''
    ... def foo():
    ...     k = 2
    ...     for i in [1, 2]:
    ...         builtins.print(k)
    ...         k = i
    ...     builtins.print(k)'''
    >>> node = ast.parse(code)
    >>> res = pm.gather(LazynessAnalysis, node)
    >>> (res['i'], res['k']) == (sys.maxsize, 1)
    True
    >>> code = '''
    ... def foo():
    ...     k = 2
    ...     for i in [1, 2]:
    ...         builtins.print(k)
    ...         k = i
    ...         builtins.print(k)'''
    >>> node = ast.parse(code)
    >>> res = pm.gather(LazynessAnalysis, node)
    >>> (res['i'], res['k']) == (sys.maxsize, 2)
    True
    >>> code = '''
    ... def foo():
    ...     d = 0
    ...     for i in [0, 1]:
    ...         for j in [0, 1]:
    ...             k = 1
    ...             d += k * 2
    ...     return d'''
    >>> node = ast.parse(code)
    >>> res = pm.gather(LazynessAnalysis, node)
    >>> res['k']
    1
    >>> code = '''
    ... def foo():
    ...     k = 2
    ...     for i in [1, 2]:
    ...         builtins.print(k)'''
    >>> node = ast.parse(code)
    >>> res = pm.gather(LazynessAnalysis, node)
    >>> res['k'] == sys.maxsize
    True
    >>> code = '''
    ... def foo():
    ...     k = builtins.sum
    ...     builtins.print(k([1, 2]))'''
    >>> node = ast.parse(code)
    >>> res = pm.gather(LazynessAnalysis, node)
    >>> res['k']
    1
    """

    INF = float('inf')
    MANY = sys.maxsize

    # map variable with maximum count of use in the programm
    ResultType = dict

    def __init__(self):
        super().__init__()
        # map variable with current count of use
        self.name_count = dict()
        # map variable to variables needed to compute it
        self.use = dict()
        # gather variables which can't be compute later. (variables used
        # to compute it have changed
        self.dead = set()
        # count use of variable before first assignation in the loop
        # {variable: (count, is_assigned)}
        self.pre_loop_count = dict()
        # prevent any form of Forward Substitution at omp frontier
        self.in_omp = set()
        self.name_to_nodes = dict()

    def modify(self, name):
        # if we modify a variable, all variables that needed it
        # to be compute are dead and its aliases too
        dead_vars = [var for var, deps in self.use.items() if name in deps]
        self.dead.update(dead_vars)
        for var in dead_vars:
            dead_aliases = [alias.id for alias in self.name_to_nodes[var]
                            if isinstance(alias, ast.Name)]
            self.dead.update(dead_aliases)

    def assign_to(self, node, from_):
        if isinstance(node, ast.Name):
            self.name_to_nodes.setdefault(node.id, set()).add(node)
        # a reassigned variable is not dead anymore
        if node.id in self.dead:
            self.dead.remove(node.id)
        # we keep the bigger possible number of use
        self.result[node.id] = max(self.result.get(node.id, 0),
                                   self.name_count.get(node.id, 0))
        # assign variable don't come from before omp pragma anymore
        self.in_omp.discard(node.id)
        # count number of use in the loop before first reassign
        pre_loop = self.pre_loop_count.setdefault(node.id, (0, True))
        if not pre_loop[1]:
            self.pre_loop_count[node.id] = (pre_loop[0], True)
        # note this variable as modified
        self.modify(node.id)
        # prepare a new variable count
        self.name_count[node.id] = 0
        self.use[node.id] = set(from_)

    def visit(self, node):
        old_omp = self.in_omp
        omp_nodes = md.get(node, openmp.OMPDirective)
        if omp_nodes:
            self.in_omp = set(self.name_count.keys())
        super(LazynessAnalysis, self).visit(node)
        if omp_nodes:
            new_nodes = set(self.name_count).difference(self.in_omp)
            for omp_node in omp_nodes:
                for n in omp_node.deps:
                    if isinstance(n, ast.Name):
                        self.result[n.id] = LazynessAnalysis.INF
            self.dead.update(new_nodes)
        self.in_omp = old_omp

    def visit_FunctionDef(self, node):
        self.ids = self.gather(Identifiers, node)
        self.generic_visit(node)

        # update result with last name_count values
        for name, val in self.name_count.items():
            old_val = self.result.get(name, 0)
            self.result[name] = max(old_val, val)

    def visit_Assign(self, node):
        md.visit(self, node)
        if node.value:
            self.visit(node.value)
            ids = self.gather(Identifiers, node.value)
        else:
            ids = set()
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        for target in targets:
            if isinstance(target, ast.Name):
                self.assign_to(target, ids)
                if node.value not in self.pure_expressions:
                    self.result[target.id] = LazynessAnalysis.INF
            elif isinstance(target, (ast.Subscript)) or isattr(target):
                # if we modify just a part of a variable, it can't be lazy
                var_name = get_variable(target)
                if isinstance(var_name, ast.Name):
                    # variable is modified so other variables that use it dies
                    self.modify(var_name.id)
                    # and this variable can't be lazy
                    self.result[var_name.id] = LazynessAnalysis.INF
            else:
                raise PythranSyntaxError("Assign to unknown node", node)

    visit_AnnAssign = visit_Assign

    def visit_AugAssign(self, node):
        md.visit(self, node)
        # augassigned variable can't be lazy
        self.visit(node.value)
        if isinstance(node.target, ast.Name):
            # variable is modified so other variables that use it dies
            self.modify(node.target.id)
            # and this variable can't be lazy
            self.result[node.target.id] = LazynessAnalysis.INF
        elif isinstance(node.target, ast.Subscript) or isattr(node.target):
            var_name = get_variable(node.target)
            # variable is modified so other variables that use it dies
            self.modify(var_name.id)
            # and this variable can't be lazy
            self.result[var_name.id] = LazynessAnalysis.INF
        else:
            raise PythranSyntaxError("AugAssign to unknown node", node)

    def visit_Name(self, node):
        if isinstance(node.ctx, ast.Load) and node.id in self.use:
            # we only care about variable local to the function

            def is_loc_var(x):
                return isinstance(x, ast.Name) and x.id in self.ids
            alias_names = [var for var in self.aliases[node]
                           if is_loc_var(var)]
            alias_names = {x.id for x in alias_names}
            alias_names.add(node.id)
            for alias in alias_names:
                if (node.id in self.dead or
                        node.id in self.in_omp):
                    self.result[alias] = LazynessAnalysis.INF
                elif alias in self.name_count:
                    self.name_count[alias] += 1
                    # init value as pre_use variable and count it
                    pre_loop = self.pre_loop_count.setdefault(alias,
                                                              (0, False))
                    if not pre_loop[1]:
                        self.pre_loop_count[alias] = (pre_loop[0] + 1, False)
                else:
                    # a variable may alias to assigned value (with a = b, 'b'
                    # alias on 'a' as modifying 'a' will modify 'b' too)
                    pass
        elif isinstance(node.ctx, ast.Param):
            self.name_count[node.id] = 0
            self.use[node.id] = set()
        elif isinstance(node.ctx, ast.Store):
            # Store is only for exception
            self.name_count[node.id] = LazynessAnalysis.INF
            self.use[node.id] = set()
        else:
            # we ignore globals
            pass

    def visit_If(self, node):
        md.visit(self, node)
        self.visit(node.test)
        old_count = dict(self.name_count)
        old_dead = set(self.dead)
        old_deps = {a: set(b) for a, b in self.use.items()}

        # wrap body in a list if we come from an ifExp
        body = node.body if isinstance(node.body, list) else [node.body]
        for stmt in body:
            self.visit(stmt)

        mid_count = self.name_count
        mid_dead = self.dead
        mid_deps = self.use

        self.name_count = old_count
        self.dead = old_dead
        self.use = old_deps

        # wrap orelse in a list if we come from an ifExp
        orelse = (node.orelse if isinstance(node.orelse, list)
                  else [node.orelse])
        for stmt in orelse:
            self.visit(stmt)

        # merge use variable
        for key in self.use:
            if key in mid_deps:
                self.use[key].update(mid_deps[key])
        for key in mid_deps:
            if key not in self.use:
                self.use[key] = set(mid_deps[key])

        # value is the worse case of both branches
        names = set(self.name_count.keys()).union(mid_count.keys())
        for name in names:
            val_body = mid_count.get(name, 0)
            val_else = self.name_count.get(name, 0)
            self.name_count[name] = max(val_body, val_else)

        # dead var are still dead
        self.dead.update(mid_dead)

    visit_IfExp = visit_If

    def visit_loop(self, body):
        # we start a new loop so we init the "at start of loop use" counter
        old_pre_count = self.pre_loop_count
        self.pre_loop_count = dict()

        # do visit body
        for stmt in body:
            self.visit(stmt)

        # variable use in loop but not assigned are no lazy
        no_assign = [n for n, (_, a) in self.pre_loop_count.items()
                     if not a]
        self.result.update(zip(no_assign,
                               [LazynessAnalysis.MANY] * len(no_assign)))
        # lazyness value is the max of previous lazyness and lazyness for one
        # iteration in the loop
        for k, v in self.pre_loop_count.items():
            loop_value = v[0] + self.name_count[k]
            self.result[k] = max(self.result.get(k, 0), loop_value)
        # variable dead at the end of the loop but use at the beginning of it
        # can't be lazy
        dead = self.dead.intersection(self.pre_loop_count)
        self.result.update(zip(dead, [LazynessAnalysis.INF] * len(dead)))
        # merge previous count of "use at start of loop" and current state.
        for k, v in old_pre_count.items():
            if v[1] or k not in self.pre_loop_count:
                self.pre_loop_count[k] = v
            else:
                self.pre_loop_count[k] = (v[0] + self.pre_loop_count[k][0],
                                          self.pre_loop_count[k][1])

    def visit_For(self, node):
        md.visit(self, node)
        ids = self.gather(Identifiers, node.iter)
        if isinstance(node.target, ast.Name):
            self.assign_to(node.target, ids)
            self.result[node.target.id] = LazynessAnalysis.INF
        else:
            err = "Assignation in for loop not to a Name"
            raise PythranSyntaxError(err, node)

        self.visit_loop(node.body)

        for stmt in node.orelse:
            self.visit(stmt)

    def visit_While(self, node):
        md.visit(self, node)
        self.visit(node.test)

        self.visit_loop(node.body)

        for stmt in node.orelse:
            self.visit(stmt)

    def func_args_lazyness(self, func_name, args, node):
        for fun in self.aliases[func_name]:
            if isinstance(fun, ast.Call):  # call to partial functions
                self.func_args_lazyness(fun.args[0], fun.args[1:] + args, node)
            elif fun in self.argument_effects:
                # when there is an argument effect, apply "modify" to the arg
                for i, arg in enumerate(self.argument_effects[fun]):
                    # check len of args as default is 11 args
                    if arg and len(args) > i:
                        if isinstance(args[i], ast.Name):
                            self.modify(args[i].id)
            elif isinstance(fun, ast.Name):
                # it may be a variable to a function. Lazyness will be compute
                # correctly thanks to aliasing
                continue
            else:
                # conservative choice
                for arg in args:
                    self.modify(arg)

    def visit_Call(self, node):
        """
        Compute use of variables in a function call.

        Each arg is use once and function name too.
        Information about modified arguments is forwarded to
        func_args_lazyness.
        """
        md.visit(self, node)
        for arg in node.args:
            self.visit(arg)
        self.func_args_lazyness(node.func, node.args, node)
        self.visit(node.func)
