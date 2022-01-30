from typing import Optional, Dict, Union
from mypy.types import TypeVarLikeType, TypeVarType, ParamSpecType, ParamSpecFlavor
from mypy.nodes import ParamSpecExpr, TypeVarExpr, TypeVarLikeExpr, SymbolTableNode


class TypeVarLikeScope:
    """Scope that holds bindings for type variables and parameter specifications.

    Node fullname -> TypeVarLikeType.
    """

    def __init__(self,
                 parent: 'Optional[TypeVarLikeScope]' = None,
                 is_class_scope: bool = False,
                 prohibited: 'Optional[TypeVarLikeScope]' = None) -> None:
        """Initializer for TypeVarLikeScope

        Parameters:
          parent: the outer scope for this scope
          is_class_scope: True if this represents a generic class
          prohibited: Type variables that aren't strictly in scope exactly,
                      but can't be bound because they're part of an outer class's scope.
        """
        self.scope: Dict[str, TypeVarLikeType] = {}
        self.parent = parent
        self.func_id = 0
        self.class_id = 0
        self.is_class_scope = is_class_scope
        self.prohibited = prohibited
        if parent is not None:
            self.func_id = parent.func_id
            self.class_id = parent.class_id

    def get_function_scope(self) -> 'Optional[TypeVarLikeScope]':
        """Get the nearest parent that's a function scope, not a class scope"""
        it: Optional[TypeVarLikeScope] = self
        while it is not None and it.is_class_scope:
            it = it.parent
        return it

    def allow_binding(self, fullname: str) -> bool:
        if fullname in self.scope:
            return False
        elif self.parent and not self.parent.allow_binding(fullname):
            return False
        elif self.prohibited and not self.prohibited.allow_binding(fullname):
            return False
        return True

    def method_frame(self) -> 'TypeVarLikeScope':
        """A new scope frame for binding a method"""
        return TypeVarLikeScope(self, False, None)

    def class_frame(self) -> 'TypeVarLikeScope':
        """A new scope frame for binding a class. Prohibits *this* class's tvars"""
        return TypeVarLikeScope(self.get_function_scope(), True, self)

    def bind_new(self, name: str, tvar_expr: TypeVarLikeExpr) -> TypeVarLikeType:
        if self.is_class_scope:
            self.class_id += 1
            i = self.class_id
        else:
            self.func_id -= 1
            i = self.func_id
        if isinstance(tvar_expr, TypeVarExpr):
            tvar_def: TypeVarLikeType = TypeVarType(
                name,
                tvar_expr.fullname,
                i,
                values=tvar_expr.values,
                upper_bound=tvar_expr.upper_bound,
                variance=tvar_expr.variance,
                line=tvar_expr.line,
                column=tvar_expr.column
            )
        elif isinstance(tvar_expr, ParamSpecExpr):
            tvar_def = ParamSpecType(
                name,
                tvar_expr.fullname,
                i,
                flavor=ParamSpecFlavor.BARE,
                upper_bound=tvar_expr.upper_bound,
                line=tvar_expr.line,
                column=tvar_expr.column
            )
        else:
            assert False
        self.scope[tvar_expr.fullname] = tvar_def
        return tvar_def

    def bind_existing(self, tvar_def: TypeVarLikeType) -> None:
        self.scope[tvar_def.fullname] = tvar_def

    def get_binding(self, item: Union[str, SymbolTableNode]) -> Optional[TypeVarLikeType]:
        fullname = item.fullname if isinstance(item, SymbolTableNode) else item
        assert fullname is not None
        if fullname in self.scope:
            return self.scope[fullname]
        elif self.parent is not None:
            return self.parent.get_binding(fullname)
        else:
            return None

    def __str__(self) -> str:
        me = ", ".join('{}: {}`{}'.format(k, v.name, v.id) for k, v in self.scope.items())
        if self.parent is None:
            return me
        return "{} <- {}".format(str(self.parent), me)
