"""Test cases for the constraint solver used in type inference."""

from __future__ import annotations

from mypy.constraints import SUBTYPE_OF, SUPERTYPE_OF, Constraint
from mypy.solve import Bounds, Graph, solve_constraints, transitive_closure
from mypy.test.helpers import Suite, assert_equal
from mypy.test.typefixture import TypeFixture
from mypy.types import Type, TypeVarId, TypeVarLikeType, TypeVarType


class SolveSuite(Suite):
    def setUp(self) -> None:
        self.fx = TypeFixture()

    def test_empty_input(self) -> None:
        self.assert_solve([], [], [])

    def test_simple_supertype_constraints(self) -> None:
        self.assert_solve([self.fx.t], [self.supc(self.fx.t, self.fx.a)], [self.fx.a])
        self.assert_solve(
            [self.fx.t],
            [self.supc(self.fx.t, self.fx.a), self.supc(self.fx.t, self.fx.b)],
            [self.fx.a],
        )

    def test_simple_subtype_constraints(self) -> None:
        self.assert_solve([self.fx.t], [self.subc(self.fx.t, self.fx.a)], [self.fx.a])
        self.assert_solve(
            [self.fx.t],
            [self.subc(self.fx.t, self.fx.a), self.subc(self.fx.t, self.fx.b)],
            [self.fx.b],
        )

    def test_both_kinds_of_constraints(self) -> None:
        self.assert_solve(
            [self.fx.t],
            [self.supc(self.fx.t, self.fx.b), self.subc(self.fx.t, self.fx.a)],
            [self.fx.b],
        )

    def test_unsatisfiable_constraints(self) -> None:
        # The constraints are impossible to satisfy.
        self.assert_solve(
            [self.fx.t],
            [self.supc(self.fx.t, self.fx.a), self.subc(self.fx.t, self.fx.b)],
            [None],
        )

    def test_exactly_specified_result(self) -> None:
        self.assert_solve(
            [self.fx.t],
            [self.supc(self.fx.t, self.fx.b), self.subc(self.fx.t, self.fx.b)],
            [self.fx.b],
        )

    def test_multiple_variables(self) -> None:
        self.assert_solve(
            [self.fx.t, self.fx.s],
            [
                self.supc(self.fx.t, self.fx.b),
                self.supc(self.fx.s, self.fx.c),
                self.subc(self.fx.t, self.fx.a),
            ],
            [self.fx.b, self.fx.c],
        )

    def test_no_constraints_for_var(self) -> None:
        self.assert_solve([self.fx.t], [], [self.fx.uninhabited])
        self.assert_solve(
            [self.fx.t, self.fx.s], [], [self.fx.uninhabited, self.fx.uninhabited]
        )
        self.assert_solve(
            [self.fx.t, self.fx.s],
            [self.supc(self.fx.s, self.fx.a)],
            [self.fx.uninhabited, self.fx.a],
        )

    def test_simple_constraints_with_dynamic_type(self) -> None:
        self.assert_solve(
            [self.fx.t], [self.supc(self.fx.t, self.fx.anyt)], [self.fx.anyt]
        )
        self.assert_solve(
            [self.fx.t],
            [self.supc(self.fx.t, self.fx.anyt), self.supc(self.fx.t, self.fx.anyt)],
            [self.fx.anyt],
        )
        self.assert_solve(
            [self.fx.t],
            [self.supc(self.fx.t, self.fx.anyt), self.supc(self.fx.t, self.fx.a)],
            [self.fx.anyt],
        )

        self.assert_solve(
            [self.fx.t], [self.subc(self.fx.t, self.fx.anyt)], [self.fx.anyt]
        )
        self.assert_solve(
            [self.fx.t],
            [self.subc(self.fx.t, self.fx.anyt), self.subc(self.fx.t, self.fx.anyt)],
            [self.fx.anyt],
        )
        # self.assert_solve([self.fx.t],
        #                   [self.subc(self.fx.t, self.fx.anyt),
        #                    self.subc(self.fx.t, self.fx.a)],
        #                   [self.fx.anyt])
        # TODO: figure out what this should be after changes to meet(any, X)

    def test_both_normal_and_any_types_in_results(self) -> None:
        # If one of the bounds is any, we promote the other bound to
        # any as well, since otherwise the type range does not make sense.
        self.assert_solve(
            [self.fx.t],
            [self.supc(self.fx.t, self.fx.a), self.subc(self.fx.t, self.fx.anyt)],
            [self.fx.anyt],
        )

        self.assert_solve(
            [self.fx.t],
            [self.supc(self.fx.t, self.fx.anyt), self.subc(self.fx.t, self.fx.a)],
            [self.fx.anyt],
        )

    def test_poly_no_constraints(self) -> None:
        self.assert_solve(
            [self.fx.t, self.fx.u],
            [],
            [self.fx.uninhabited, self.fx.uninhabited],
            allow_polymorphic=True,
        )

    def test_poly_trivial_free(self) -> None:
        self.assert_solve(
            [self.fx.t, self.fx.u],
            [self.subc(self.fx.t, self.fx.a)],
            [self.fx.a, self.fx.u],
            [self.fx.u],
            allow_polymorphic=True,
        )

    def test_poly_free_pair(self) -> None:
        self.assert_solve(
            [self.fx.t, self.fx.u],
            [self.subc(self.fx.t, self.fx.u)],
            [self.fx.t, self.fx.t],
            [self.fx.t],
            allow_polymorphic=True,
        )

    def test_poly_free_pair_with_bounds(self) -> None:
        t_prime = self.fx.t.copy_modified(upper_bound=self.fx.b)
        self.assert_solve(
            [self.fx.t, self.fx.ub],
            [self.subc(self.fx.t, self.fx.ub)],
            [t_prime, t_prime],
            [t_prime],
            allow_polymorphic=True,
        )

    def test_poly_free_pair_with_bounds_uninhabited(self) -> None:
        self.assert_solve(
            [self.fx.ub, self.fx.uc],
            [self.subc(self.fx.ub, self.fx.uc)],
            [self.fx.uninhabited, self.fx.uninhabited],
            [],
            allow_polymorphic=True,
        )

    def test_poly_bounded_chain(self) -> None:
        # B <: T <: U <: S <: A
        self.assert_solve(
            [self.fx.t, self.fx.u, self.fx.s],
            [
                self.supc(self.fx.t, self.fx.b),
                self.subc(self.fx.t, self.fx.u),
                self.subc(self.fx.u, self.fx.s),
                self.subc(self.fx.s, self.fx.a),
            ],
            [self.fx.b, self.fx.b, self.fx.b],
            allow_polymorphic=True,
        )

    def test_poly_reverse_overlapping_chain(self) -> None:
        # A :> T <: S :> B
        self.assert_solve(
            [self.fx.t, self.fx.s],
            [
                self.subc(self.fx.t, self.fx.s),
                self.subc(self.fx.t, self.fx.a),
                self.supc(self.fx.s, self.fx.b),
            ],
            [self.fx.a, self.fx.a],
            allow_polymorphic=True,
        )

    def test_poly_reverse_split_chain(self) -> None:
        # B :> T <: S :> A
        self.assert_solve(
            [self.fx.t, self.fx.s],
            [
                self.subc(self.fx.t, self.fx.s),
                self.subc(self.fx.t, self.fx.b),
                self.supc(self.fx.s, self.fx.a),
            ],
            [self.fx.b, self.fx.a],
            allow_polymorphic=True,
        )

    def test_poly_unsolvable_chain(self) -> None:
        # A <: T <: U <: S <: B
        self.assert_solve(
            [self.fx.t, self.fx.u, self.fx.s],
            [
                self.supc(self.fx.t, self.fx.a),
                self.subc(self.fx.t, self.fx.u),
                self.subc(self.fx.u, self.fx.s),
                self.subc(self.fx.s, self.fx.b),
            ],
            [None, None, None],
            allow_polymorphic=True,
        )

    def test_simple_chain_closure(self) -> None:
        self.assert_transitive_closure(
            [self.fx.t.id, self.fx.s.id],
            [
                self.supc(self.fx.t, self.fx.b),
                self.subc(self.fx.t, self.fx.s),
                self.subc(self.fx.s, self.fx.a),
            ],
            {(self.fx.t.id, self.fx.s.id)},
            {self.fx.t.id: {self.fx.b}, self.fx.s.id: {self.fx.b}},
            {self.fx.t.id: {self.fx.a}, self.fx.s.id: {self.fx.a}},
        )

    def test_reverse_chain_closure(self) -> None:
        self.assert_transitive_closure(
            [self.fx.t.id, self.fx.s.id],
            [
                self.subc(self.fx.t, self.fx.s),
                self.subc(self.fx.t, self.fx.a),
                self.supc(self.fx.s, self.fx.b),
            ],
            {(self.fx.t.id, self.fx.s.id)},
            {self.fx.t.id: set(), self.fx.s.id: {self.fx.b}},
            {self.fx.t.id: {self.fx.a}, self.fx.s.id: set()},
        )

    def test_secondary_constraint_closure(self) -> None:
        self.assert_transitive_closure(
            [self.fx.t.id, self.fx.s.id],
            [self.supc(self.fx.s, self.fx.gt), self.subc(self.fx.s, self.fx.ga)],
            set(),
            {self.fx.t.id: set(), self.fx.s.id: {self.fx.gt}},
            {self.fx.t.id: {self.fx.a}, self.fx.s.id: {self.fx.ga}},
        )

    def assert_solve(
        self,
        vars: list[TypeVarLikeType],
        constraints: list[Constraint],
        results: list[None | Type],
        free_vars: list[TypeVarLikeType] | None = None,
        allow_polymorphic: bool = False,
    ) -> None:
        if free_vars is None:
            free_vars = []
        actual, actual_free = solve_constraints(
            vars, constraints, allow_polymorphic=allow_polymorphic
        )
        assert_equal(actual, results)
        assert_equal(actual_free, free_vars)

    def assert_transitive_closure(
        self,
        vars: list[TypeVarId],
        constraints: list[Constraint],
        graph: Graph,
        lowers: Bounds,
        uppers: Bounds,
    ) -> None:
        actual_graph, actual_lowers, actual_uppers = transitive_closure(
            vars, constraints
        )
        # Add trivial elements.
        for v in vars:
            graph.add((v, v))
        assert_equal(actual_graph, graph)
        assert_equal(dict(actual_lowers), lowers)
        assert_equal(dict(actual_uppers), uppers)

    def supc(self, type_var: TypeVarType, bound: Type) -> Constraint:
        return Constraint(type_var, SUPERTYPE_OF, bound)

    def subc(self, type_var: TypeVarType, bound: Type) -> Constraint:
        return Constraint(type_var, SUBTYPE_OF, bound)
