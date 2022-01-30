"""Test cases for the constraint solver used in type inference."""

from typing import List, Union, Tuple, Optional

from mypy.test.helpers import Suite, assert_equal
from mypy.constraints import SUPERTYPE_OF, SUBTYPE_OF, Constraint
from mypy.solve import solve_constraints
from mypy.test.typefixture import TypeFixture
from mypy.types import Type, TypeVarType, TypeVarId


class SolveSuite(Suite):
    def setUp(self) -> None:
        self.fx = TypeFixture()

    def test_empty_input(self) -> None:
        self.assert_solve([], [], [])

    def test_simple_supertype_constraints(self) -> None:
        self.assert_solve([self.fx.t.id],
                          [self.supc(self.fx.t, self.fx.a)],
                          [(self.fx.a, self.fx.o)])
        self.assert_solve([self.fx.t.id],
                          [self.supc(self.fx.t, self.fx.a),
                           self.supc(self.fx.t, self.fx.b)],
                          [(self.fx.a, self.fx.o)])

    def test_simple_subtype_constraints(self) -> None:
        self.assert_solve([self.fx.t.id],
                          [self.subc(self.fx.t, self.fx.a)],
                          [self.fx.a])
        self.assert_solve([self.fx.t.id],
                          [self.subc(self.fx.t, self.fx.a),
                           self.subc(self.fx.t, self.fx.b)],
                          [self.fx.b])

    def test_both_kinds_of_constraints(self) -> None:
        self.assert_solve([self.fx.t.id],
                          [self.supc(self.fx.t, self.fx.b),
                           self.subc(self.fx.t, self.fx.a)],
                          [(self.fx.b, self.fx.a)])

    def test_unsatisfiable_constraints(self) -> None:
        # The constraints are impossible to satisfy.
        self.assert_solve([self.fx.t.id],
                          [self.supc(self.fx.t, self.fx.a),
                           self.subc(self.fx.t, self.fx.b)],
                          [None])

    def test_exactly_specified_result(self) -> None:
        self.assert_solve([self.fx.t.id],
                          [self.supc(self.fx.t, self.fx.b),
                           self.subc(self.fx.t, self.fx.b)],
                          [(self.fx.b, self.fx.b)])

    def test_multiple_variables(self) -> None:
        self.assert_solve([self.fx.t.id, self.fx.s.id],
                          [self.supc(self.fx.t, self.fx.b),
                           self.supc(self.fx.s, self.fx.c),
                           self.subc(self.fx.t, self.fx.a)],
                          [(self.fx.b, self.fx.a), (self.fx.c, self.fx.o)])

    def test_no_constraints_for_var(self) -> None:
        self.assert_solve([self.fx.t.id],
                          [],
                          [self.fx.uninhabited])
        self.assert_solve([self.fx.t.id, self.fx.s.id],
                          [],
                          [self.fx.uninhabited, self.fx.uninhabited])
        self.assert_solve([self.fx.t.id, self.fx.s.id],
                          [self.supc(self.fx.s, self.fx.a)],
                          [self.fx.uninhabited, (self.fx.a, self.fx.o)])

    def test_simple_constraints_with_dynamic_type(self) -> None:
        self.assert_solve([self.fx.t.id],
                          [self.supc(self.fx.t, self.fx.anyt)],
                          [(self.fx.anyt, self.fx.anyt)])
        self.assert_solve([self.fx.t.id],
                          [self.supc(self.fx.t, self.fx.anyt),
                           self.supc(self.fx.t, self.fx.anyt)],
                          [(self.fx.anyt, self.fx.anyt)])
        self.assert_solve([self.fx.t.id],
                          [self.supc(self.fx.t, self.fx.anyt),
                           self.supc(self.fx.t, self.fx.a)],
                          [(self.fx.anyt, self.fx.anyt)])

        self.assert_solve([self.fx.t.id],
                          [self.subc(self.fx.t, self.fx.anyt)],
                          [(self.fx.anyt, self.fx.anyt)])
        self.assert_solve([self.fx.t.id],
                          [self.subc(self.fx.t, self.fx.anyt),
                           self.subc(self.fx.t, self.fx.anyt)],
                          [(self.fx.anyt, self.fx.anyt)])
        # self.assert_solve([self.fx.t.id],
        #                   [self.subc(self.fx.t, self.fx.anyt),
        #                    self.subc(self.fx.t, self.fx.a)],
        #                   [(self.fx.anyt, self.fx.anyt)])
        # TODO: figure out what this should be after changes to meet(any, X)

    def test_both_normal_and_any_types_in_results(self) -> None:
        # If one of the bounds is any, we promote the other bound to
        # any as well, since otherwise the type range does not make sense.
        self.assert_solve([self.fx.t.id],
                          [self.supc(self.fx.t, self.fx.a),
                           self.subc(self.fx.t, self.fx.anyt)],
                          [(self.fx.anyt, self.fx.anyt)])

        self.assert_solve([self.fx.t.id],
                          [self.supc(self.fx.t, self.fx.anyt),
                           self.subc(self.fx.t, self.fx.a)],
                          [(self.fx.anyt, self.fx.anyt)])

    def assert_solve(self,
                     vars: List[TypeVarId],
                     constraints: List[Constraint],
                     results: List[Union[None, Type, Tuple[Type, Type]]],
                     ) -> None:
        res: List[Optional[Type]] = []
        for r in results:
            if isinstance(r, tuple):
                res.append(r[0])
            else:
                res.append(r)
        actual = solve_constraints(vars, constraints)
        assert_equal(str(actual), str(res))

    def supc(self, type_var: TypeVarType, bound: Type) -> Constraint:
        return Constraint(type_var.id, SUPERTYPE_OF, bound)

    def subc(self, type_var: TypeVarType, bound: Type) -> Constraint:
        return Constraint(type_var.id, SUBTYPE_OF, bound)
