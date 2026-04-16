""" ParallelMaps detects parallel map(...).  """

from pythran.analyses.aliases import Aliases
from pythran.analyses.pure_expressions import PureExpressions
from pythran.passmanager import ModuleAnalysis
from pythran.tables import MODULES


class ParallelMaps(ModuleAnalysis[PureExpressions, Aliases]):

    """Yields the est of maps that could be parallel."""

    ResultType = set

    def visit_Call(self, node):
        if all(alias == MODULES['builtins']['map']
               for alias in self.aliases[node.func]):
            if all(f in self.pure_expressions
                   for f in self.aliases[node.args[0]]):
                self.result.add(node)

    def display(self, data):
        for node in data:
            print("I:", "{0} {1}".format(
                "call to the `map' intrinsic could be parallel",
                "(line {0})".format(node.lineno)
                ))
