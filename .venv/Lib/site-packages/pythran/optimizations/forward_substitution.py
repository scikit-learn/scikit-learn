"""
Replace variable that can be lazy evaluated and used only once by their full
computation code.
"""

from pythran.analyses import LazynessAnalysis, UseDefChains, DefUseChains
from pythran.analyses import Literals, Ancestors, Identifiers, CFG, IsAssigned
from pythran.passmanager import Transformation
import pythran.graph as graph

from collections import defaultdict
from copy import deepcopy
import gast as ast

try:
    from math import isfinite
except ImportError:
    from math import isinf, isnan

    def isfinite(x):
        return not isinf(x) and not isnan(x)


class Remover(ast.NodeTransformer):

    def __init__(self, nodes):
        self.nodes = nodes

    def visit_Assign(self, node):
        if node in self.nodes:
            to_prune = self.nodes[node]
            node.targets = [tgt for tgt in node.targets if tgt not in to_prune]
            if node.targets:
                return node
            else:
                return ast.Pass()
        return node

    def visit_AnnAssign(self, node):
        if node in self.nodes:
            to_prune = self.nodes[node]
            if node.target in to_prune:
                return ast.Pass()
            else:
                return node
        return node


class ForwardSubstitution(Transformation[LazynessAnalysis, UseDefChains, DefUseChains,
                                         Ancestors, CFG, Literals]):

    """
    Replace variable that can be computed later.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> pm = passmanager.PassManager("test")

    >>> node = ast.parse("def foo(): a = [2, 3]; builtins.print(a)")
    >>> _, node = pm.apply(ForwardSubstitution, node)
    >>> print(pm.dump(backend.Python, node))
    def foo():
        pass
        builtins.print([2, 3])

    >>> node = ast.parse("def foo(): a = 2; builtins.print(a + a)")
    >>> _, node = pm.apply(ForwardSubstitution, node)
    >>> print(pm.dump(backend.Python, node))
    def foo():
        a = 2
        builtins.print((2 + 2))

    >>> node = ast.parse("def foo():\\n a=b=2\\n while a: a -= 1\\n return b")
    >>> _, node = pm.apply(ForwardSubstitution, node)
    >>> print(pm.dump(backend.Python, node))
    def foo():
        a = 2
        while a:
            a -= 1
        return 2
    """

    def __init__(self):
        """ Satisfy dependencies on others analyses. """
        super().__init__()
        self.to_remove = None

    def visit_FunctionDef(self, node):
        self.to_remove = defaultdict(list)
        self.locals = self.def_use_chains.locals[node]

        # prune some assignment as a second phase, as an assignment could be
        # forward-substituted several times (in the case of constants)
        self.generic_visit(node)
        Remover(self.to_remove).visit(node)
        return node

    def visit_AugAssign(self, node):
        # Separate case for augassign, where we only handle situation where the
        # def is a literal (any kind) with a single use (this AugAssign).
        self.generic_visit(node)
        if not isinstance(node.target, ast.Name):
            return node

        defs = self.use_def_chains[node.target]
        name_defs = [dnode for dnode in defs if isinstance(dnode.node, ast.Name)]
        if len(name_defs) != 1:
            return node
        name_def, = name_defs
        dnode = name_def.node
        parent = self.ancestors[dnode][-1]
        if not isinstance(parent, ast.Assign):
            return node

        try:
            ast.literal_eval(parent.value)
        except ValueError:
            return node

        self.update = True
        return ast.Assign([node.target], ast.BinOp(deepcopy(parent.value), node.op, node.value))

    def visit_Name(self, node):
        if not isinstance(node.ctx, ast.Load):
            return node

        # OpenMP metdata are not handled by beniget, which is fine in our case
        if node not in self.use_def_chains:
            if __debug__:
                from pythran.openmp import OMPDirective
                assert any(isinstance(p, OMPDirective)
                           for p in self.ancestors[node])
            return node
        defuses = self.use_def_chains[node]

        if len(defuses) != 1:
            return node

        defuse = defuses[0]

        dnode = defuse.node
        if not isinstance(dnode, ast.Name):
            return node

        # multiple definition, which one should we forward?
        if sum(1 for d in self.locals if d.name() == dnode.id) > 1:
            return node

        # either a constant or a value
        fwd = (dnode in self.literals and
               isfinite(self.lazyness_analysis[dnode.id]))
        fwd |= self.lazyness_analysis[dnode.id] == 1

        if not fwd:
            return node

        parent = self.ancestors[dnode][-1]
        if isinstance(parent, ast.Assign):
            value = parent.value
            if dnode in self.literals:
                self.update = True
                if len(defuse.users()) == 1:
                    self.to_remove[parent].append(dnode)
                    return value
                else:
                    return value
            elif len(parent.targets) == 1:
                ids = self.gather(Identifiers, value)
                node_stmt = next(s for s in self.ancestors[node][::-1]
                                 if isinstance(s, ast.stmt))
                # Check if there is a path from `parent' to `node_stmt' that
                # modifies any of the identifier from `value'. If so, cancel the
                # forward substitution.
                worklist = [node_stmt]
                visited = {parent}
                while worklist:
                    workitem = worklist.pop()
                    if workitem in visited:
                        continue
                    visited.add(workitem)
                    for pred in self.cfg.predecessors(workitem):
                        if not graph.has_path(self.cfg, parent, pred):
                            continue

                        assigned_ids = {n.id
                                        for n in self.gather(IsAssigned,
                                                             pred)}
                        if not ids.isdisjoint(assigned_ids):
                            return node  # cancel
                        worklist.append(pred)
                self.update = True
                self.to_remove[parent].append(dnode)
                return value

        return node

class PreInliningForwardSubstitution(ForwardSubstitution):

    """
    Replace variable that can be computed later, but only if this leads to a
    one-liner that's going to be a great inlining candidate.
    """

    def visit_FunctionDef(self, node):
        # Only handle trivial cases, because this can lead to more inlining
        # opportunities.
        if all(isinstance(s, (ast.Return, ast.Assign)) for s in node.body):
            r = super(PreInliningForwardSubstitution,
                         self).visit_FunctionDef(node)
            return r
        return node
