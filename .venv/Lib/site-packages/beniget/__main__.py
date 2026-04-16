import sys

import gast as ast
from beniget import Ancestors, DefUseChains

class Beniget(ast.NodeVisitor):
    def __init__(self, filename, module):
        super(Beniget, self).__init__()

        self.filename = filename or "<stdin>"

        self.ancestors = Ancestors()
        self.ancestors.visit(module)

        self.defuses = DefUseChains(self.filename)
        self.defuses.visit(module)

        self.visit(module)

    def check_unused(self, node, skipped_types=()):
        for local_def in self.defuses.locals[node]:
            if not local_def.users():
                if local_def.name() == "_":
                    continue  # typical naming by-pass
                if isinstance(local_def.node, skipped_types):
                    continue

                location = local_def.node
                while not hasattr(location, "lineno"):
                    location = self.ancestors.parent(location)

                if isinstance(location, ast.ImportFrom):
                    if location.module == "__future__":
                        continue

                print(
                    "W: '{}' is defined but not used at {}:{}:{}".format(
                        local_def.name(),
                        self.filename,
                        location.lineno,
                        location.col_offset,
                    )
                )

    def visit_Module(self, node):
        self.generic_visit(node)
        if self.filename.endswith("__init__.py"):
            return
        self.check_unused(
            node, skipped_types=(ast.FunctionDef, ast.AsyncFunctionDef,
                                 ast.ClassDef, ast.Name)
        )

    def visit_FunctionDef(self, node):
        self.generic_visit(node)
        self.check_unused(node)

paths = sys.argv[1:] or (None,)

for path in paths:
    with open(path) if path else sys.stdin as target:
        module = ast.parse(target.read())
        Beniget(path, module)

