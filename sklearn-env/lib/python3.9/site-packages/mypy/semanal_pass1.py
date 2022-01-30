"""Block/import reachability analysis."""

from mypy.nodes import (
    MypyFile, AssertStmt, IfStmt, Block, AssignmentStmt, ExpressionStmt, ReturnStmt, ForStmt,
    Import, ImportAll, ImportFrom, ClassDef, FuncDef
)
from mypy.traverser import TraverserVisitor
from mypy.options import Options
from mypy.reachability import infer_reachability_of_if_statement, assert_will_always_fail


class SemanticAnalyzerPreAnalysis(TraverserVisitor):
    """Analyze reachability of blocks and imports and other local things.

    This runs before semantic analysis, so names have not been bound. Imports are
    also not resolved yet, so we can only access the current module.

    This determines static reachability of blocks and imports due to version and
    platform checks, among others.

    The main entry point is 'visit_file'.

    Reachability of imports needs to be determined very early in the build since
    this affects which modules will ultimately be processed.

    Consider this example:

      import sys

      def do_stuff():
          # type: () -> None:
          if sys.python_version < (3,):
              import xyz  # Only available in Python 2
              xyz.whatever()
          ...

    The block containing 'import xyz' is unreachable in Python 3 mode. The import
    shouldn't be processed in Python 3 mode, even if the module happens to exist.
    """

    def visit_file(self, file: MypyFile, fnam: str, mod_id: str, options: Options) -> None:
        self.pyversion = options.python_version
        self.platform = options.platform
        self.cur_mod_id = mod_id
        self.cur_mod_node = file
        self.options = options
        self.is_global_scope = True

        for i, defn in enumerate(file.defs):
            defn.accept(self)
            if isinstance(defn, AssertStmt) and assert_will_always_fail(defn, options):
                # We've encountered an assert that's always false,
                # e.g. assert sys.platform == 'lol'.  Truncate the
                # list of statements.  This mutates file.defs too.
                del file.defs[i + 1:]
                break

    def visit_func_def(self, node: FuncDef) -> None:
        old_global_scope = self.is_global_scope
        self.is_global_scope = False
        super().visit_func_def(node)
        self.is_global_scope = old_global_scope
        file_node = self.cur_mod_node
        if (self.is_global_scope
                and file_node.is_stub
                and node.name == '__getattr__'
                and file_node.is_package_init_file()):
            # __init__.pyi with __getattr__ means that any submodules are assumed
            # to exist, even if there is no stub. Note that we can't verify that the
            # return type is compatible, since we haven't bound types yet.
            file_node.is_partial_stub_package = True

    def visit_class_def(self, node: ClassDef) -> None:
        old_global_scope = self.is_global_scope
        self.is_global_scope = False
        super().visit_class_def(node)
        self.is_global_scope = old_global_scope

    def visit_import_from(self, node: ImportFrom) -> None:
        node.is_top_level = self.is_global_scope
        super().visit_import_from(node)

    def visit_import_all(self, node: ImportAll) -> None:
        node.is_top_level = self.is_global_scope
        super().visit_import_all(node)

    def visit_import(self, node: Import) -> None:
        node.is_top_level = self.is_global_scope
        super().visit_import(node)

    def visit_if_stmt(self, s: IfStmt) -> None:
        infer_reachability_of_if_statement(s, self.options)
        for expr in s.expr:
            expr.accept(self)
        for node in s.body:
            node.accept(self)
        if s.else_body:
            s.else_body.accept(self)

    def visit_block(self, b: Block) -> None:
        if b.is_unreachable:
            return
        super().visit_block(b)

    # The remaining methods are an optimization: don't visit nested expressions
    # of common statements, since they can have no effect.

    def visit_assignment_stmt(self, s: AssignmentStmt) -> None:
        pass

    def visit_expression_stmt(self, s: ExpressionStmt) -> None:
        pass

    def visit_return_stmt(self, s: ReturnStmt) -> None:
        pass

    def visit_for_stmt(self, s: ForStmt) -> None:
        s.body.accept(self)
        if s.else_body is not None:
            s.else_body.accept(self)
