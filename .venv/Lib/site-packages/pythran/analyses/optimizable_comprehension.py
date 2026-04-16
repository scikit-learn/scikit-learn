"""
OptimizableComp finds whether a comprehension can be optimized.
"""

from pythran.analyses.identifiers import Identifiers
from pythran.passmanager import NodeAnalysis


class OptimizableComprehension(NodeAnalysis[Identifiers]):
    """Find whether a comprehension can be optimized."""
    ResultType = set

    def check_comprehension(self, iters):
        targets = {gen.target.id for gen in iters}
        optimizable = True

        for it in iters:
            ids = self.gather(Identifiers, it)
            optimizable &= all(((ident == it.target.id) |
                                (ident not in targets)) for ident in ids)

        return optimizable

    def visit_ListComp(self, node):
        if (self.check_comprehension(node.generators)):
            self.result.add(node)

    def visit_GeneratorExp(self, node):
        if (self.check_comprehension(node.generators)):
            self.result.add(node)
