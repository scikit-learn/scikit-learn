""" FalsePolymorphism try to rename variable to avoid false polymorphism."""

from pythran.passmanager import Transformation
from pythran.analyses import DefUseChains, UseDefChains, Identifiers

import gast as ast
import re


class FalsePolymorphism(Transformation[DefUseChains, UseDefChains]):

    """
    Rename variable when possible to avoid false polymorphism.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("def foo(): a = 12; a = 'babar'")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(FalsePolymorphism, node)
    >>> print(pm.dump(backend.Python, node))
    def foo():
        a = 12
        a_ = 'babar'
    """

    def visit_FunctionDef(self, node):

        # reset available identifier names
        # removing local identifiers from the list so that first occurrence can
        # actually use the slot
        identifiers = self.gather(Identifiers, node)
        captured_identifiers = set()
        captured_identifiers_pattern = re.compile('^__pythran_boxed_(?:args_)?(.*)$')
        for def_ in self.def_use_chains.locals[node]:
            match = captured_identifiers_pattern.match(def_.name())
            if match:
                captured_identifiers.add(match.group(1))
            try:
                identifiers.remove(def_.name())
            except KeyError:
                pass

        # compute all reachable nodes from each def. This builds a bag of def
        # that should have the same name
        visited_defs = set()
        for def_ in self.def_use_chains.locals[node]:
            if def_.name() in captured_identifiers:
                continue
            if def_ in visited_defs:
                continue

            associated_defs = set()

            # fill the bag of associated defs, going through users and defs
            to_process = [def_]
            while to_process:
                curr = to_process.pop()
                if curr in associated_defs:
                    continue
                if curr.name() != def_.name():
                    continue
                associated_defs.add(curr)
                for u in curr.users():
                    to_process.append(u)
                curr_udc = (d for d in self.use_def_chains.get(curr.node, [])
                            if isinstance(d.node, ast.Name))
                to_process.extend(curr_udc)

            visited_defs.update(associated_defs)

            # find a new identifier
            local_identifier = def_.name()
            name = local_identifier
            while name in identifiers:
                name += "_"
            identifiers.add(name)

            # don't rename first candidate
            if name == local_identifier:
                continue

            # actual renaming of each node in the bag
            self.update = True
            for d in associated_defs:
                dn = d.node
                if isinstance(dn, ast.Name) and dn.id == local_identifier:
                    dn.id = name
        return node
