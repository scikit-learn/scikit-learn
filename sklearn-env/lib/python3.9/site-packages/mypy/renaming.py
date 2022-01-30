from contextlib import contextmanager
from typing import Dict, Iterator, List
from typing_extensions import Final

from mypy.nodes import (
    Block, AssignmentStmt, NameExpr, MypyFile, FuncDef, Lvalue, ListExpr, TupleExpr,
    WhileStmt, ForStmt, BreakStmt, ContinueStmt, TryStmt, WithStmt, StarExpr, ImportFrom,
    MemberExpr, IndexExpr, Import, ClassDef
)
from mypy.traverser import TraverserVisitor

# Scope kinds
FILE: Final = 0
FUNCTION: Final = 1
CLASS: Final = 2


class VariableRenameVisitor(TraverserVisitor):
    """Rename variables to allow redefinition of variables.

    For example, consider this code:

      x = 0
      f(x)

      x = "a"
      g(x)

    It will be transformed like this:

      x' = 0
      f(x')

      x = "a"
      g(x)

    There will be two independent variables (x' and x) that will have separate
    inferred types. The publicly exposed variant will get the non-suffixed name.
    This is the last definition at module top level and the first definition
    (argument) within a function.

    Renaming only happens for assignments within the same block. Renaming is
    performed before semantic analysis, immediately after parsing.

    The implementation performs a rudimentary static analysis. The analysis is
    overly conservative to keep things simple.
    """

    def __init__(self) -> None:
        # Counter for labeling new blocks
        self.block_id = 0
        # Number of surrounding try statements that disallow variable redefinition
        self.disallow_redef_depth = 0
        # Number of surrounding loop statements
        self.loop_depth = 0
        # Map block id to loop depth.
        self.block_loop_depth: Dict[int, int] = {}
        # Stack of block ids being processed.
        self.blocks: List[int] = []
        # List of scopes; each scope maps short (unqualified) name to block id.
        self.var_blocks: List[Dict[str, int]] = []

        # References to variables that we may need to rename. List of
        # scopes; each scope is a mapping from name to list of collections
        # of names that refer to the same logical variable.
        self.refs: List[Dict[str, List[List[NameExpr]]]] = []
        # Number of reads of the most recent definition of a variable (per scope)
        self.num_reads: List[Dict[str, int]] = []
        # Kinds of nested scopes (FILE, FUNCTION or CLASS)
        self.scope_kinds: List[int] = []

    def visit_mypy_file(self, file_node: MypyFile) -> None:
        """Rename variables within a file.

        This is the main entry point to this class.
        """
        self.clear()
        with self.enter_scope(FILE), self.enter_block():
            for d in file_node.defs:
                d.accept(self)

    def visit_func_def(self, fdef: FuncDef) -> None:
        # Conservatively do not allow variable defined before a function to
        # be redefined later, since function could refer to either definition.
        self.reject_redefinition_of_vars_in_scope()

        with self.enter_scope(FUNCTION), self.enter_block():
            for arg in fdef.arguments:
                name = arg.variable.name
                # 'self' can't be redefined since it's special as it allows definition of
                # attributes. 'cls' can't be used to define attributes so we can ignore it.
                can_be_redefined = name != 'self'  # TODO: Proper check
                self.record_assignment(arg.variable.name, can_be_redefined)
                self.handle_arg(name)

            for stmt in fdef.body.body:
                stmt.accept(self)

    def visit_class_def(self, cdef: ClassDef) -> None:
        self.reject_redefinition_of_vars_in_scope()
        with self.enter_scope(CLASS):
            super().visit_class_def(cdef)

    def visit_block(self, block: Block) -> None:
        with self.enter_block():
            super().visit_block(block)

    def visit_while_stmt(self, stmt: WhileStmt) -> None:
        with self.enter_loop():
            super().visit_while_stmt(stmt)

    def visit_for_stmt(self, stmt: ForStmt) -> None:
        stmt.expr.accept(self)
        self.analyze_lvalue(stmt.index, True)
        # Also analyze as non-lvalue so that every for loop index variable is assumed to be read.
        stmt.index.accept(self)
        with self.enter_loop():
            stmt.body.accept(self)
        if stmt.else_body:
            stmt.else_body.accept(self)

    def visit_break_stmt(self, stmt: BreakStmt) -> None:
        self.reject_redefinition_of_vars_in_loop()

    def visit_continue_stmt(self, stmt: ContinueStmt) -> None:
        self.reject_redefinition_of_vars_in_loop()

    def visit_try_stmt(self, stmt: TryStmt) -> None:
        # Variables defined by a try statement get special treatment in the
        # type checker which allows them to be always redefined, so no need to
        # do renaming here.
        with self.enter_try():
            super().visit_try_stmt(stmt)

    def visit_with_stmt(self, stmt: WithStmt) -> None:
        for expr in stmt.expr:
            expr.accept(self)
        for target in stmt.target:
            if target is not None:
                self.analyze_lvalue(target)
        # We allow redefinitions in the body of a with statement for
        # convenience.  This is unsafe since with statements can affect control
        # flow by catching exceptions, but this is rare except for
        # assertRaises() and other similar functions, where the exception is
        # raised by the last statement in the body, which usually isn't a
        # problem.
        stmt.body.accept(self)

    def visit_import(self, imp: Import) -> None:
        for id, as_id in imp.ids:
            self.record_assignment(as_id or id, False)

    def visit_import_from(self, imp: ImportFrom) -> None:
        for id, as_id in imp.names:
            self.record_assignment(as_id or id, False)

    def visit_assignment_stmt(self, s: AssignmentStmt) -> None:
        s.rvalue.accept(self)
        for lvalue in s.lvalues:
            self.analyze_lvalue(lvalue)

    def analyze_lvalue(self, lvalue: Lvalue, is_nested: bool = False) -> None:
        """Process assignment; in particular, keep track of (re)defined names.

        Args:
            is_nested: True for non-outermost Lvalue in a multiple assignment such as
                "x, y = ..."
        """
        if isinstance(lvalue, NameExpr):
            name = lvalue.name
            is_new = self.record_assignment(name, True)
            if is_new:
                self.handle_def(lvalue)
            else:
                self.handle_refine(lvalue)
            if is_nested:
                # This allows these to be redefined freely even if never read. Multiple
                # assignment like "x, _ _ = y" defines dummy variables that are never read.
                self.handle_ref(lvalue)
        elif isinstance(lvalue, (ListExpr, TupleExpr)):
            for item in lvalue.items:
                self.analyze_lvalue(item, is_nested=True)
        elif isinstance(lvalue, MemberExpr):
            lvalue.expr.accept(self)
        elif isinstance(lvalue, IndexExpr):
            lvalue.base.accept(self)
            lvalue.index.accept(self)
        elif isinstance(lvalue, StarExpr):
            # Propagate is_nested since in a typical use case like "x, *rest = ..." 'rest' may
            # be freely reused.
            self.analyze_lvalue(lvalue.expr, is_nested=is_nested)

    def visit_name_expr(self, expr: NameExpr) -> None:
        self.handle_ref(expr)

    # Helpers for renaming references

    def handle_arg(self, name: str) -> None:
        """Store function argument."""
        self.refs[-1][name] = [[]]
        self.num_reads[-1][name] = 0

    def handle_def(self, expr: NameExpr) -> None:
        """Store new name definition."""
        name = expr.name
        names = self.refs[-1].setdefault(name, [])
        names.append([expr])
        self.num_reads[-1][name] = 0

    def handle_refine(self, expr: NameExpr) -> None:
        """Store assignment to an existing name (that replaces previous value, if any)."""
        name = expr.name
        if name in self.refs[-1]:
            names = self.refs[-1][name]
            if not names:
                names.append([])
            names[-1].append(expr)

    def handle_ref(self, expr: NameExpr) -> None:
        """Store reference to defined name."""
        name = expr.name
        if name in self.refs[-1]:
            names = self.refs[-1][name]
            if not names:
                names.append([])
            names[-1].append(expr)
        num_reads = self.num_reads[-1]
        num_reads[name] = num_reads.get(name, 0) + 1

    def flush_refs(self) -> None:
        """Rename all references within the current scope.

        This will be called at the end of a scope.
        """
        is_func = self.scope_kinds[-1] == FUNCTION
        for name, refs in self.refs[-1].items():
            if len(refs) == 1:
                # Only one definition -- no renaming needed.
                continue
            if is_func:
                # In a function, don't rename the first definition, as it
                # may be an argument that must preserve the name.
                to_rename = refs[1:]
            else:
                # At module top level, don't rename the final definition,
                # as it will be publicly visible outside the module.
                to_rename = refs[:-1]
            for i, item in enumerate(to_rename):
                self.rename_refs(item, i)
        self.refs.pop()

    def rename_refs(self, names: List[NameExpr], index: int) -> None:
        name = names[0].name
        new_name = name + "'" * (index + 1)
        for expr in names:
            expr.name = new_name

    # Helpers for determining which assignments define new variables

    def clear(self) -> None:
        self.blocks = []
        self.var_blocks = []

    @contextmanager
    def enter_block(self) -> Iterator[None]:
        self.block_id += 1
        self.blocks.append(self.block_id)
        self.block_loop_depth[self.block_id] = self.loop_depth
        try:
            yield
        finally:
            self.blocks.pop()

    @contextmanager
    def enter_try(self) -> Iterator[None]:
        self.disallow_redef_depth += 1
        try:
            yield
        finally:
            self.disallow_redef_depth -= 1

    @contextmanager
    def enter_loop(self) -> Iterator[None]:
        self.loop_depth += 1
        try:
            yield
        finally:
            self.loop_depth -= 1

    def current_block(self) -> int:
        return self.blocks[-1]

    @contextmanager
    def enter_scope(self, kind: int) -> Iterator[None]:
        self.var_blocks.append({})
        self.refs.append({})
        self.num_reads.append({})
        self.scope_kinds.append(kind)
        try:
            yield
        finally:
            self.flush_refs()
            self.var_blocks.pop()
            self.num_reads.pop()
            self.scope_kinds.pop()

    def is_nested(self) -> int:
        return len(self.var_blocks) > 1

    def reject_redefinition_of_vars_in_scope(self) -> None:
        """Make it impossible to redefine defined variables in the current scope.

        This is used if we encounter a function definition that
        can make it ambiguous which definition is live. Example:

          x = 0

          def f() -> int:
              return x

          x = ''  # Error -- cannot redefine x across function definition
        """
        var_blocks = self.var_blocks[-1]
        for key in var_blocks:
            var_blocks[key] = -1

    def reject_redefinition_of_vars_in_loop(self) -> None:
        """Reject redefinition of variables in the innermost loop.

        If there is an early exit from a loop, there may be ambiguity about which
        value may escape the loop. Example where this matters:

          while f():
              x = 0
              if g():
                  break
              x = ''  # Error -- not a redefinition
          reveal_type(x)  # int

        This method ensures that the second assignment to 'x' doesn't introduce a new
        variable.
        """
        var_blocks = self.var_blocks[-1]
        for key, block in var_blocks.items():
            if self.block_loop_depth.get(block) == self.loop_depth:
                var_blocks[key] = -1

    def record_assignment(self, name: str, can_be_redefined: bool) -> bool:
        """Record assignment to given name and return True if it defines a new variable.

        Args:
            can_be_redefined: If True, allows assignment in the same block to redefine
                this name (if this is a new definition)
        """
        if self.num_reads[-1].get(name, -1) == 0:
            # Only set, not read, so no reason to redefine
            return False
        if self.disallow_redef_depth > 0:
            # Can't redefine within try/with a block.
            can_be_redefined = False
        block = self.current_block()
        var_blocks = self.var_blocks[-1]
        if name not in var_blocks:
            # New definition in this scope.
            if can_be_redefined:
                # Store the block where this was defined to allow redefinition in
                # the same block only.
                var_blocks[name] = block
            else:
                # This doesn't support arbitrary redefinition.
                var_blocks[name] = -1
            return True
        elif var_blocks[name] == block:
            # Redefinition -- defines a new variable with the same name.
            return True
        else:
            # Assigns to an existing variable.
            return False
