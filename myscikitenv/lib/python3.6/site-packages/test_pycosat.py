import sys
import copy
import random
from os.path import basename
import unittest

import pycosat
from pycosat import solve, itersolve

# -------------------------- utility functions ---------------------------

def read_cnf(path):
    """
    read a DIMACS cnf formatted file from `path`, and return the clauses
    and number of variables
    """
    clauses = []
    for line in open(path):
        parts = line.split()
        if not parts or parts[0] == 'c':
            continue
        if parts[0] == 'p':
            assert len(parts) == 4
            assert parts[1] == 'cnf'
            n_vars, n_clauses = [int(n) for n in parts[2:4]]
            continue
        if parts[0] == '%':
            break
        assert parts[-1] == '0'
        clauses.append([int(lit) for lit in parts[:-1]])
    assert len(clauses) == n_clauses
    return clauses, n_vars

def evaluate(clauses, sol):
    """
    evaluate the clauses with the solution
    """
    sol_vars = {} # variable number -> bool
    for i in sol:
        sol_vars[abs(i)] = bool(i > 0)
    return all(any(sol_vars[abs(i)] ^ bool(i < 0) for i in clause)
               for clause in clauses)

def py_itersolve(clauses):
    while True:
        sol = pycosat.solve(clauses)
        if isinstance(sol, list):
            yield sol
            clauses.append([-x for x in sol])
        else: # no more solutions -- stop iteration
            return

def process_cnf_file(path):
    sys.stdout.write('%30s:  ' % basename(path))
    sys.stdout.flush()

    clauses, n_vars = read_cnf(path)
    sys.stdout.write('vars: %6d   cls: %6d   ' % (n_vars, len(clauses)))
    sys.stdout.flush()
    n_sol = 0
    for sol in itersolve(clauses, n_vars):
        sys.stdout.write('.')
        sys.stdout.flush()
        assert evaluate(clauses, sol)
        n_sol += 1
    sys.stdout.write("%d\n" % n_sol)
    sys.stdout.flush()
    return n_sol

# -------------------------- test clauses --------------------------------

# p cnf 5 3
# 1 -5 4 0
# -1 5 3 4 0
# -3 -4 0
nvars1, clauses1 = 5, [[1, -5, 4], [-1, 5, 3, 4], [-3, -4]]

# p cnf 2 2
# -1 0
# 1 0
nvars2, clauses2 = 2, [[-1], [1]]

# p cnf 2 3
# -1 2 0
# -1 -2 0
# 1 -2 0
nvars3, clauses3 = 2, [[-1, 2], [-1, -2], [1, -2]]

# -------------------------- actual unit tests ---------------------------

tests = []

class TestSolve(unittest.TestCase):

    def test_wrong_args(self):
        self.assertRaises(TypeError, solve, [[1, 2], [-3]], 'A')
        self.assertRaises(TypeError, solve, 1)
        self.assertRaises(TypeError, solve, 1.0)
        self.assertRaises(TypeError, solve, object())
        self.assertRaises(TypeError, solve, ['a'])
        self.assertRaises(TypeError, solve, [[1, 2], [3, None]], 5)
        self.assertRaises(ValueError, solve, [[1, 2], [3, 0]])

    def test_no_clauses(self):
        for n in range(7):
            self.assertEqual(solve([], n), [-i for i in range(1, n + 1)])

    def test_cnf1(self):
        self.assertEqual(solve(clauses1), [1, -2, -3, -4, 5])
        if sys.version_info[0] == 2:
            cls = [[long(lit) for lit in clause] for clause in clauses1]
            self.assertEqual(solve(cls), [1, -2, -3, -4, 5])

    def test_iter_clauses(self):
        self.assertEqual(solve(iter(clauses1)), [1, -2, -3, -4, 5])

    def test_each_clause_iter(self):
        self.assertEqual(solve([iter(clause) for clause in clauses1]),
                         [1, -2, -3, -4, 5])

    def test_tuple_caluses(self):
        self.assertEqual(solve(tuple(clauses1)), [1, -2, -3, -4, 5])

    def test_each_clause_tuples(self):
        self.assertEqual(solve([tuple(clause) for clause in clauses1]),
                         [1, -2, -3, -4, 5])

    def test_gen_clauses(self):
        def gen_clauses():
            for clause in clauses1:
                yield clause
        self.assertEqual(solve(gen_clauses()), [1, -2, -3, -4, 5])

    def test_each_clause_gen(self):
        self.assertEqual(solve([(x for x in clause) for clause in clauses1]),
                         [1, -2, -3, -4, 5])

    def test_bad_iter(self):
        class Liar:
            def __iter__(self): return None
        self.assertRaises(TypeError, solve, Liar())

    def test_cnf2(self):
        self.assertEqual(solve(clauses2), "UNSAT")

    def test_cnf3(self):
        self.assertEqual(solve(clauses3), [-1, -2])

    def test_cnf3_3vars(self):
        self.assertEqual(solve(clauses3, vars=3), [-1, -2, -3])

    def test_cnf1_prop_limit(self):
        for lim in range(1, 20):
            self.assertEqual(solve(clauses1, prop_limit=lim),
                             "UNKNOWN" if lim < 8 else [1, -2, -3, -4, 5])

    def test_cnf1_vars(self):
        self.assertEqual(solve(clauses1, vars=7),
                         [1, -2, -3, -4, 5, -6, -7])

tests.append(TestSolve)

# -----

class TestIterSolve(unittest.TestCase):

    def test_wrong_args(self):
        self.assertRaises(TypeError, itersolve, [[1, 2], [-3]], 'A')
        self.assertRaises(TypeError, itersolve, 1)
        self.assertRaises(TypeError, itersolve, 1.0)
        self.assertRaises(TypeError, itersolve, object())
        self.assertRaises(TypeError, itersolve, ['a'])
        self.assertRaises(TypeError, itersolve, [[1, 2], [3, None]], 5)
        self.assertRaises(ValueError, itersolve, [[1, 2], [3, 0]])

    def test_no_clauses(self):
        for n in range(7):
            self.assertEqual(len(list(itersolve([], vars=n))), 2 ** n)

    def test_iter_clauses(self):
        self.assertTrue(all(evaluate(clauses1, sol) for sol in
                            itersolve(iter(clauses1))))

    def test_each_clause_iter(self):
        self.assertTrue(all(evaluate(clauses1, sol) for sol in
                            itersolve([iter(clause) for clause in clauses1])))

    def test_tuple_caluses(self):
        self.assertTrue(all(evaluate(clauses1, sol) for sol in
                            itersolve(tuple(clauses1))))

    def test_each_clause_tuples(self):
        self.assertTrue(all(evaluate(clauses1, sol) for sol in
                            itersolve([tuple(clause) for clause in clauses1])))

    def test_gen_clauses(self):
        def gen_clauses():
            for clause in clauses1:
                yield clause
        self.assertTrue(all(evaluate(clauses1, sol) for sol in
                            itersolve(gen_clauses())))

    def test_each_clause_gen(self):
        self.assertTrue(all(evaluate(clauses1, sol) for sol in
                            itersolve([(x for x in clause) for clause in
                                       clauses1])))

    def test_bad_iter(self):
        class Liar:
            def __iter__(self): return None
        self.assertRaises(TypeError, itersolve, Liar())

    def test_cnf1(self):
        for sol in itersolve(clauses1, nvars1):
            #sys.stderr.write('%r\n' % repr(sol))
            self.assertTrue(evaluate(clauses1, sol))

        sols = list(itersolve(clauses1, vars=nvars1))
        self.assertEqual(len(sols), 18)
        # ensure solutions are unique
        self.assertEqual(len(set(tuple(sol) for sol in sols)), 18)

    def test_shuffle_clauses(self):
        ref_sols = set(tuple(sol) for sol in itersolve(clauses1))
        for _ in range(10):
            cnf = copy.deepcopy(clauses1)
            # shuffling the clauses does not change the solutions
            random.shuffle(cnf)
            self.assertEqual(set(tuple(sol) for sol in itersolve(cnf)),
                             ref_sols)

    def test_many_clauses(self):
        ref_sols = set(tuple(sol) for sol in itersolve(clauses1))
        # repeating the clauses many times does not change the solutions
        cnf = 100 * copy.deepcopy(clauses1)
        self.assertEqual(set(tuple(sol) for sol in itersolve(cnf)),
                         ref_sols)

    def test_cnf2(self):
        self.assertEqual(list(itersolve(clauses2, nvars2)), [])

    def test_cnf3_3vars(self):
        self.assertEqual(list(itersolve(clauses3, 3)),
                         [[-1, -2, -3], [-1, -2, 3]])

    def test_cnf1_prop_limit(self):
        self.assertEqual(list(itersolve(clauses1, prop_limit=2)), [])

tests.append(TestIterSolve)

# ------------------------------------------------------------------------

def run(verbosity=1, repeat=1):
    print("sys.prefix: %s" % sys.prefix)
    print("sys.version: %s" % sys.version)
    print("pycosat version: %r" % pycosat.__version__)

    suite = unittest.TestSuite()
    for cls in tests:
        for _ in range(repeat):
            suite.addTest(unittest.makeSuite(cls))

    runner = unittest.TextTestRunner(verbosity=verbosity)
    return runner.run(suite)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        run()
    else:
        for path in sys.argv[1:]:
            process_cnf_file(path)
