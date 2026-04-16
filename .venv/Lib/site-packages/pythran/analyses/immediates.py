"""
Immediates gathers immediates. For now, only integers within shape are
and argument of functions flagged as immediate_arguments are
considered as immediates
"""

from pythran.analyses import Aliases
from pythran.passmanager import NodeAnalysis
from pythran.utils import pythran_builtin, isnum, ispowi

_make_shape = pythran_builtin('make_shape')


class Immediates(NodeAnalysis[Aliases]):
    ResultType = set

    def visit_BinOp(self, node):
        self.generic_visit(node)
        if ispowi(node):
            self.result.add(node.right)

    def visit_AugAssign(self, node):
        self.generic_visit(node)
        if ispowi(node):
            self.result.add(node.value)

    def visit_Call(self, node):
        func_aliases = self.aliases[node.func]
        for alias in func_aliases:
            if getattr(alias, "immediate_arguments", []):
                for i, arg in enumerate(node.args):
                    if i in alias.immediate_arguments:
                        self.result.add(arg)

        if len(func_aliases) == 1 and next(iter(func_aliases)) is _make_shape:
            self.result.update(a for a in node.args
                               if isnum(a)
                               and isinstance(a.value, int)
                               and a.value >= 0)
            return

        return self.generic_visit(node)
