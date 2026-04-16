"""
Scope computes scope information
"""

from pythran.analyses.ancestors import AncestorsWithBody
from pythran.analyses.use_def_chain import DefUseChains
from pythran.passmanager import FunctionAnalysis

from collections import defaultdict
import gast as ast


def all_users(d):
    """
    Gather users of a definition, including users of an augmented assign.
    """
    visited = set()
    dname = d.name()
    def all_users_impl(d):
        if d in visited:
            return
        visited.add(d)
        for u in d.users():
            if u.name() == dname:
                yield u
                yield from all_users_impl(u)
    return all_users_impl(d)


class Scope(FunctionAnalysis[AncestorsWithBody, DefUseChains]):
    '''
    Associate each variable declaration with the node that defines it

    Whenever possible, associate the variable declaration to an assignment,
    otherwise to a node that defines a bloc (e.g. a For)
    This takes OpenMP information into accounts!
    The result is a dictionary with nodes as key and set of names as values
    '''

    ResultType = lambda: defaultdict(set)
    def __init__(self):
        super().__init__()
        self.decl_holders = (ast.FunctionDef, ast.For,
                             ast.excepthandler,
                             ast.While, ast.If, tuple)

    def visit_OMPDirective(self, node):
        for dep in node.deps:
            if dep in node.private_deps:
                continue
            if isinstance(dep, ast.Name):
                self.openmp_deps.setdefault(dep.id, []).append(dep)

    def visit_FunctionDef(self, node):
        self.ancestors = self.ancestors_with_body

        # first gather some info about OpenMP declarations
        self.openmp_deps = dict()
        self.generic_visit(node)

        name_to_defs = dict()
        for def_ in self.def_use_chains.locals[node]:
            name_to_defs.setdefault(def_.name(), []).append(def_)

        # then compute scope informations
        # unlike use-def chains, this takes OpenMP annotations into account
        for name, defs in name_to_defs.items():
            # get all refs to that name
            refs = [d.node for d in defs] + [u.node
                                             for d in defs for u in all_users(d)]
            # add OpenMP refs (well, the parent of the holding stmt)
            refs.extend(self.ancestors[d][-3]   # -3 to get the right parent
                        for d in self.openmp_deps.get(name, []))
            # get their ancestors
            ancestors = [self.ancestors[ref] for ref in refs]
            # common ancestors
            prefixes = [p for p in zip(*ancestors) if len(set(p)) == 1]
            common = prefixes[-1][0]  # the last common ancestor

            # now try to attach the scope to an assignment.
            # This will be the first assignment found in the bloc
            if isinstance(common, self.decl_holders):
                # get all refs that define that name
                refs = [d.node for d in defs]
                refs.extend(self.openmp_deps.get(name, []))

                # get their parent
                prefs = set()
                for r in refs:
                    ancestor = r
                    # walk up the ancestor tree until we find the one
                    # right before common
                    while self.ancestors[ancestor][-1] is not common:
                        ancestor = self.ancestors[ancestor][-1]
                    prefs.add(ancestor)

                # set the defining statement to the first assign in the body
                # unless another statements uses it before
                # or the common itselfs holds a dependency
                if common not in prefs:
                    body = common if isinstance(common, tuple) else common.body
                    for c in body:
                        if c in prefs:
                            if isinstance(c, ast.Assign):
                                common = c
                            break
            self.result[common].add(name)
