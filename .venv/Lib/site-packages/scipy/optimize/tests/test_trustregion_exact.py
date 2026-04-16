"""
Unit tests for trust-region iterative subproblem.

"""
import pytest
import numpy as np
from scipy.optimize._trustregion_exact import (
    estimate_smallest_singular_value,
    singular_leading_submatrix,
    IterativeSubproblem)
from scipy.linalg import (svd, get_lapack_funcs, det, qr, norm)
from numpy.testing import (assert_array_equal,
                           assert_equal, assert_array_almost_equal)


def random_entry(n, min_eig, max_eig, case, rng=None):
    rng = np.random.default_rng(rng)

    # Generate random matrix
    rand = rng.uniform(low=-1, high=1, size=(n, n))

    # QR decomposition
    Q, _, _ = qr(rand, pivoting='True')

    # Generate random eigenvalues
    eigvalues = rng.uniform(low=min_eig, high=max_eig, size=n)
    eigvalues = np.sort(eigvalues)[::-1]

    # Generate matrix
    Qaux = np.multiply(eigvalues, Q)
    A = np.dot(Qaux, Q.T)

    # Generate gradient vector accordingly
    # to the case is being tested.
    if case == 'hard':
        g = np.zeros(n)
        g[:-1] = rng.uniform(low=-1, high=1, size=n-1)
        g = np.dot(Q, g)
    elif case == 'jac_equal_zero':
        g = np.zeros(n)
    else:
        g = rng.uniform(low=-1, high=1, size=n)

    return A, g


class TestEstimateSmallestSingularValue:

    def test_for_ill_condiotioned_matrix(self):

        # Ill-conditioned triangular matrix
        C = np.array([[1, 2, 3, 4],
                      [0, 0.05, 60, 7],
                      [0, 0, 0.8, 9],
                      [0, 0, 0, 10]])

        # Get svd decomposition
        U, s, Vt = svd(C)

        # Get smallest singular value and correspondent right singular vector.
        smin_svd = s[-1]
        zmin_svd = Vt[-1, :]

        # Estimate smallest singular value
        smin, zmin = estimate_smallest_singular_value(C)

        # Check the estimation
        assert_array_almost_equal(smin, smin_svd, decimal=8)
        assert_array_almost_equal(abs(zmin), abs(zmin_svd), decimal=8)


class TestSingularLeadingSubmatrix:

    def test_for_already_singular_leading_submatrix(self):

        # Define test matrix A.
        # Note that the leading 2x2 submatrix is singular.
        A = np.array([[1, 2, 3],
                      [2, 4, 5],
                      [3, 5, 6]])

        # Get Cholesky from lapack functions
        cholesky, = get_lapack_funcs(('potrf',), (A,))

        # Compute Cholesky Decomposition
        c, k = cholesky(A, lower=False, overwrite_a=False, clean=True)

        delta, v = singular_leading_submatrix(A, c, k)

        A[k-1, k-1] += delta

        # Check if the leading submatrix is singular.
        assert_array_almost_equal(det(A[:k, :k]), 0)

        # Check if `v` fulfil the specified properties
        quadratic_term = np.dot(v, np.dot(A, v))
        assert_array_almost_equal(quadratic_term, 0)

    def test_for_simetric_indefinite_matrix(self):

        # Define test matrix A.
        # Note that the leading 5x5 submatrix is indefinite.
        A = np.asarray([[1, 2, 3, 7, 8],
                        [2, 5, 5, 9, 0],
                        [3, 5, 11, 1, 2],
                        [7, 9, 1, 7, 5],
                        [8, 0, 2, 5, 8]])

        # Get Cholesky from lapack functions
        cholesky, = get_lapack_funcs(('potrf',), (A,))

        # Compute Cholesky Decomposition
        c, k = cholesky(A, lower=False, overwrite_a=False, clean=True)

        delta, v = singular_leading_submatrix(A, c, k)

        A[k-1, k-1] += delta

        # Check if the leading submatrix is singular.
        assert_array_almost_equal(det(A[:k, :k]), 0)

        # Check if `v` fulfil the specified properties
        quadratic_term = np.dot(v, np.dot(A, v))
        assert_array_almost_equal(quadratic_term, 0)

    def test_for_first_element_equal_to_zero(self):

        # Define test matrix A.
        # Note that the leading 2x2 submatrix is singular.
        A = np.array([[0, 3, 11],
                      [3, 12, 5],
                      [11, 5, 6]])

        # Get Cholesky from lapack functions
        cholesky, = get_lapack_funcs(('potrf',), (A,))

        # Compute Cholesky Decomposition
        c, k = cholesky(A, lower=False, overwrite_a=False, clean=True)

        delta, v = singular_leading_submatrix(A, c, k)

        A[k-1, k-1] += delta

        # Check if the leading submatrix is singular
        assert_array_almost_equal(det(A[:k, :k]), 0)

        # Check if `v` fulfil the specified properties
        quadratic_term = np.dot(v, np.dot(A, v))
        assert_array_almost_equal(quadratic_term, 0)


class TestIterativeSubproblem:

    def test_for_the_easy_case(self):

        # `H` is chosen such that `g` is not orthogonal to the
        # eigenvector associated with the smallest eigenvalue `s`.
        H = [[10, 2, 3, 4],
             [2, 1, 7, 1],
             [3, 7, 1, 7],
             [4, 1, 7, 2]]
        g = [1, 1, 1, 1]

        # Trust Radius
        trust_radius = 1

        # Solve Subproblem
        subprob = IterativeSubproblem(x=0,
                                      fun=lambda x: 0,
                                      jac=lambda x: np.array(g),
                                      hess=lambda x: np.array(H),
                                      k_easy=1e-10,
                                      k_hard=1e-10)
        p, hits_boundary = subprob.solve(trust_radius)

        assert_array_almost_equal(p, [0.00393332, -0.55260862,
                                      0.67065477, -0.49480341])
        assert_array_almost_equal(hits_boundary, True)

    def test_for_the_hard_case(self):

        # `H` is chosen such that `g` is orthogonal to the
        # eigenvector associated with the smallest eigenvalue `s`.
        H = [[10, 2, 3, 4],
             [2, 1, 7, 1],
             [3, 7, 1, 7],
             [4, 1, 7, 2]]
        g = [6.4852641521327437, 1, 1, 1]
        s = -8.2151519874416614

        # Trust Radius
        trust_radius = 1

        # Solve Subproblem
        subprob = IterativeSubproblem(x=0,
                                      fun=lambda x: 0,
                                      jac=lambda x: np.array(g),
                                      hess=lambda x: np.array(H),
                                      k_easy=1e-10,
                                      k_hard=1e-10)
        p, hits_boundary = subprob.solve(trust_radius)

        assert_array_almost_equal(-s, subprob.lambda_current)

    def test_for_interior_convergence(self):

        H = [[1.812159, 0.82687265, 0.21838879, -0.52487006, 0.25436988],
             [0.82687265, 2.66380283, 0.31508988, -0.40144163, 0.08811588],
             [0.21838879, 0.31508988, 2.38020726, -0.3166346, 0.27363867],
             [-0.52487006, -0.40144163, -0.3166346, 1.61927182, -0.42140166],
             [0.25436988, 0.08811588, 0.27363867, -0.42140166, 1.33243101]]

        g = [0.75798952, 0.01421945, 0.33847612, 0.83725004, -0.47909534]

        # Solve Subproblem
        subprob = IterativeSubproblem(x=0,
                                      fun=lambda x: 0,
                                      jac=lambda x: np.array(g),
                                      hess=lambda x: np.array(H))
        p, hits_boundary = subprob.solve(1.1)

        assert_array_almost_equal(p, [-0.68585435, 0.1222621, -0.22090999,
                                      -0.67005053, 0.31586769])
        assert_array_almost_equal(hits_boundary, False)
        assert_array_almost_equal(subprob.lambda_current, 0)
        assert_array_almost_equal(subprob.niter, 1)

    def test_for_jac_equal_zero(self):

        H = [[0.88547534, 2.90692271, 0.98440885, -0.78911503, -0.28035809],
             [2.90692271, -0.04618819, 0.32867263, -0.83737945, 0.17116396],
             [0.98440885, 0.32867263, -0.87355957, -0.06521957, -1.43030957],
             [-0.78911503, -0.83737945, -0.06521957, -1.645709, -0.33887298],
             [-0.28035809, 0.17116396, -1.43030957, -0.33887298, -1.68586978]]

        g = [0, 0, 0, 0, 0]

        # Solve Subproblem
        subprob = IterativeSubproblem(x=0,
                                      fun=lambda x: 0,
                                      jac=lambda x: np.array(g),
                                      hess=lambda x: np.array(H),
                                      k_easy=1e-10,
                                      k_hard=1e-10)
        p, hits_boundary = subprob.solve(1.1)

        assert_array_almost_equal(p, [0.06910534, -0.01432721,
                                      -0.65311947, -0.23815972,
                                      -0.84954934])
        assert_array_almost_equal(hits_boundary, True)

    def test_for_jac_very_close_to_zero(self):

        H = [[0.88547534, 2.90692271, 0.98440885, -0.78911503, -0.28035809],
             [2.90692271, -0.04618819, 0.32867263, -0.83737945, 0.17116396],
             [0.98440885, 0.32867263, -0.87355957, -0.06521957, -1.43030957],
             [-0.78911503, -0.83737945, -0.06521957, -1.645709, -0.33887298],
             [-0.28035809, 0.17116396, -1.43030957, -0.33887298, -1.68586978]]

        g = [0, 0, 0, 0, 1e-15]

        # Solve Subproblem
        subprob = IterativeSubproblem(x=0,
                                      fun=lambda x: 0,
                                      jac=lambda x: np.array(g),
                                      hess=lambda x: np.array(H),
                                      k_easy=1e-10,
                                      k_hard=1e-10)
        p, hits_boundary = subprob.solve(1.1)

        assert_array_almost_equal(p, [0.06910534, -0.01432721,
                                      -0.65311947, -0.23815972,
                                      -0.84954934])
        assert_array_almost_equal(hits_boundary, True)

    @pytest.mark.thread_unsafe(reason="fails in parallel")
    @pytest.mark.fail_slow(10)
    def test_for_random_entries(self):
        rng = np.random.default_rng(1)

        # Dimension
        n = 5

        for case in ('easy', 'hard', 'jac_equal_zero'):

            eig_limits = [(-20, -15),
                          (-10, -5),
                          (-10, 0),
                          (-5, 5),
                          (-10, 10),
                          (0, 10),
                          (5, 10),
                          (15, 20)]

            for min_eig, max_eig in eig_limits:
                # Generate random symmetric matrix H with
                # eigenvalues between min_eig and max_eig.
                H, g = random_entry(n, min_eig, max_eig, case, rng=rng)

                # Trust radius
                trust_radius_list = [0.1, 0.3, 0.6, 0.8, 1, 1.2, 3.3, 5.5, 10]

                for trust_radius in trust_radius_list:
                    # Solve subproblem with very high accuracy
                    subprob_ac = IterativeSubproblem(0,
                                                     lambda x: 0,
                                                     lambda x: g,
                                                     lambda x: H,
                                                     k_easy=1e-10,
                                                     k_hard=1e-10)

                    p_ac, hits_boundary_ac = subprob_ac.solve(trust_radius)

                    # Compute objective function value
                    J_ac = 1/2*np.dot(p_ac, np.dot(H, p_ac))+np.dot(g, p_ac)

                    stop_criteria = [(0.1, 2),
                                     (0.5, 1.1),
                                     (0.9, 1.01)]

                    for k_opt, k_trf in stop_criteria:

                        # k_easy and k_hard computed in function
                        # of k_opt and k_trf accordingly to
                        # Conn, A. R., Gould, N. I., & Toint, P. L. (2000).
                        # "Trust region methods". Siam. p. 197.
                        k_easy = min(k_trf-1,
                                     1-np.sqrt(k_opt))
                        k_hard = 1-k_opt

                        # Solve subproblem
                        subprob = IterativeSubproblem(0,
                                                      lambda x: 0,
                                                      lambda x: g,
                                                      lambda x: H,
                                                      k_easy=k_easy,
                                                      k_hard=k_hard)
                        p, hits_boundary = subprob.solve(trust_radius)

                        # Compute objective function value
                        J = 1/2*np.dot(p, np.dot(H, p))+np.dot(g, p)

                        # Check if it respect k_trf
                        if hits_boundary:
                            assert_array_equal(np.abs(norm(p)-trust_radius) <=
                                               (k_trf-1)*trust_radius, True)
                        else:
                            assert_equal(norm(p) <= trust_radius, True)

                        # Check if it respect k_opt
                        assert_equal(J <= k_opt*J_ac, True)


    def test_for_finite_number_of_iterations(self):
        """Regression test for gh-12513"""
        H = np.array(
            [[3.67335930e01, -2.52334820e02, 1.15477558e01, -1.19933725e-03,
              -2.06408851e03, -2.05821411e00, -2.52334820e02, -6.52076924e02,
              -2.71362566e-01, -1.98885126e00, 1.22085415e00, 2.30220713e00,
              -9.71278532e-02, -5.11210123e-01, -1.00399562e00, 1.43319679e-01,
              6.03815471e00, -6.38719934e-02, 1.65623929e-01],
             [-2.52334820e02, 1.76757312e03, -9.92814996e01, 1.06533600e-02,
              1.44442941e04, 1.43811694e01, 1.76757312e03, 4.56694461e03,
              2.22263363e00, 1.62977318e01, -7.81539315e00, -1.24938012e01,
              6.74029088e-01, 3.22802671e00, 5.14978971e00, -9.58561209e-01,
              -3.92199895e01, 4.47201278e-01, -1.17866744e00],
             [1.15477558e01, -9.92814996e01, 3.63872363e03, -4.40007197e-01,
              -9.55435081e02, -1.13985105e00, -9.92814996e01, -2.58307255e02,
              -5.21335218e01, -3.77485107e02, -6.75338369e01, -1.89457169e02,
              5.67828623e00, 5.82402681e00, 1.72734354e01, -4.29114840e00,
              -7.84885258e01, 3.17594634e00, 2.45242852e00],
             [-1.19933725e-03, 1.06533600e-02, -4.40007197e-01, 5.73576663e-05,
              1.01563710e-01, 1.18838745e-04, 1.06533600e-02, 2.76535767e-02,
              6.25788669e-03, 4.50699620e-02, 8.64152333e-03, 2.27772377e-02,
              -8.51026855e-04, 1.65316383e-04, 1.38977551e-03, 5.51629259e-04,
              1.38447755e-02, -5.17956723e-04, -1.29260347e-04],
             [-2.06408851e03, 1.44442941e04, -9.55435081e02, 1.01563710e-01,
              1.23101825e05, 1.26467259e02, 1.44442941e04, 3.74590279e04,
              2.18498571e01, 1.60254460e02, -7.52977260e01, -1.17989623e02,
              6.58253160e00, 3.14949206e01, 4.98527190e01, -9.33338661e00,
              -3.80465752e02, 4.33872213e00, -1.14768816e01],
             [-2.05821411e00, 1.43811694e01, -1.13985105e00, 1.18838745e-04,
              1.26467259e02, 1.46226198e-01, 1.43811694e01, 3.74509252e01,
              2.76928748e-02, 2.03023837e-01, -8.84279903e-02, -1.29523344e-01,
              8.06424434e-03, 3.83330661e-02, 5.81579023e-02, -1.12874980e-02,
              -4.48118297e-01, 5.15022284e-03, -1.41501894e-02],
             [-2.52334820e02, 1.76757312e03, -9.92814996e01, 1.06533600e-02,
              1.44442941e04, 1.43811694e01, 1.76757312e03, 4.56694461e03,
              2.22263363e00, 1.62977318e01, -7.81539315e00, -1.24938012e01,
              6.74029088e-01, 3.22802671e00, 5.14978971e00, -9.58561209e-01,
              -3.92199895e01, 4.47201278e-01, -1.17866744e00],
             [-6.52076924e02, 4.56694461e03, -2.58307255e02, 2.76535767e-02,
              3.74590279e04, 3.74509252e01, 4.56694461e03, 1.18278398e04,
              5.82242837e00, 4.26867612e01, -2.03167952e01, -3.22894255e01,
              1.75705078e00, 8.37153730e00, 1.32246076e01, -2.49238529e00,
              -1.01316422e02, 1.16165466e00, -3.09390862e00],
             [-2.71362566e-01, 2.22263363e00, -5.21335218e01, 6.25788669e-03,
              2.18498571e01, 2.76928748e-02, 2.22263363e00, 5.82242837e00,
              4.36278066e01, 3.14836583e02, -2.04747938e01, -3.05535101e01,
              -1.24881456e-01, 1.15775394e01, 4.06907410e01, -1.39317748e00,
              -3.90902798e01, -9.71716488e-02, 1.06851340e-01],
             [-1.98885126e00, 1.62977318e01, -3.77485107e02, 4.50699620e-02,
              1.60254460e02, 2.03023837e-01, 1.62977318e01, 4.26867612e01,
              3.14836583e02, 2.27255216e03, -1.47029712e02, -2.19649109e02,
              -8.83963155e-01, 8.28571708e01, 2.91399776e02, -9.97382920e00,
              -2.81069124e02, -6.94946614e-01, 7.38151960e-01],
             [1.22085415e00, -7.81539315e00, -6.75338369e01, 8.64152333e-03,
              -7.52977260e01, -8.84279903e-02, -7.81539315e00, -2.03167952e01,
              -2.04747938e01, -1.47029712e02, 7.83372613e01, 1.64416651e02,
              -4.30243758e00, -2.59579610e01, -6.25644064e01, 6.69974667e00,
              2.31011701e02, -2.68540084e00, 5.44531151e00],
             [2.30220713e00, -1.24938012e01, -1.89457169e02, 2.27772377e-02,
              -1.17989623e02, -1.29523344e-01, -1.24938012e01, -3.22894255e01,
              -3.05535101e01, -2.19649109e02, 1.64416651e02, 3.75893031e02,
              -7.42084715e00, -4.56437599e01, -1.11071032e02, 1.18761368e01,
              4.78724142e02, -5.06804139e00, 8.81448081e00],
             [-9.71278532e-02, 6.74029088e-01, 5.67828623e00, -8.51026855e-04,
              6.58253160e00, 8.06424434e-03, 6.74029088e-01, 1.75705078e00,
              -1.24881456e-01, -8.83963155e-01, -4.30243758e00, -7.42084715e00,
              9.62009425e-01, 1.53836355e00, 2.23939458e00, -8.01872920e-01,
              -1.92191084e01, 3.77713908e-01, -8.32946970e-01],
             [-5.11210123e-01, 3.22802671e00, 5.82402681e00, 1.65316383e-04,
              3.14949206e01, 3.83330661e-02, 3.22802671e00, 8.37153730e00,
              1.15775394e01, 8.28571708e01, -2.59579610e01, -4.56437599e01,
              1.53836355e00, 2.63851056e01, 7.34859767e01, -4.39975402e00,
              -1.12015747e02, 5.11542219e-01, -2.64962727e00],
             [-1.00399562e00, 5.14978971e00, 1.72734354e01, 1.38977551e-03,
              4.98527190e01, 5.81579023e-02, 5.14978971e00, 1.32246076e01,
              4.06907410e01, 2.91399776e02, -6.25644064e01, -1.11071032e02,
              2.23939458e00, 7.34859767e01, 2.36535458e02, -1.09636675e01,
              -2.72152068e02, 6.65888059e-01, -6.29295273e00],
             [1.43319679e-01, -9.58561209e-01, -4.29114840e00, 5.51629259e-04,
              -9.33338661e00, -1.12874980e-02, -9.58561209e-01, -2.49238529e00,
              -1.39317748e00, -9.97382920e00, 6.69974667e00, 1.18761368e01,
              -8.01872920e-01, -4.39975402e00, -1.09636675e01, 1.16820748e00,
              3.00817252e01, -4.51359819e-01, 9.82625204e-01],
             [6.03815471e00, -3.92199895e01, -7.84885258e01, 1.38447755e-02,
              -3.80465752e02, -4.48118297e-01, -3.92199895e01, -1.01316422e02,
              -3.90902798e01, -2.81069124e02, 2.31011701e02, 4.78724142e02,
              -1.92191084e01, -1.12015747e02, -2.72152068e02, 3.00817252e01,
              1.13232557e03, -1.33695932e01, 2.22934659e01],
             [-6.38719934e-02, 4.47201278e-01, 3.17594634e00, -5.17956723e-04,
              4.33872213e00, 5.15022284e-03, 4.47201278e-01, 1.16165466e00,
              -9.71716488e-02, -6.94946614e-01, -2.68540084e00, -5.06804139e00,
              3.77713908e-01, 5.11542219e-01, 6.65888059e-01, -4.51359819e-01,
              -1.33695932e01, 4.27994168e-01, -5.09020820e-01],
             [1.65623929e-01, -1.17866744e00, 2.45242852e00, -1.29260347e-04,
              -1.14768816e01, -1.41501894e-02, -1.17866744e00, -3.09390862e00,
              1.06851340e-01, 7.38151960e-01, 5.44531151e00, 8.81448081e00,
              -8.32946970e-01, -2.64962727e00, -6.29295273e00, 9.82625204e-01,
              2.22934659e01, -5.09020820e-01, 4.09964606e00]]
        )
        J = np.array([
            -2.53298102e-07, 1.76392040e-06, 1.74776130e-06, -4.19479903e-10,
            1.44167498e-05, 1.41703911e-08, 1.76392030e-06, 4.96030153e-06,
            -2.35771675e-07, -1.68844985e-06, 4.29218258e-07, 6.65445159e-07,
            -3.87045830e-08, -3.17236594e-07, -1.21120169e-06, 4.59717313e-08,
            1.67123246e-06, 1.46624675e-08, 4.22723383e-08
        ])

        subproblem_maxiter = [None, 10]
        for maxiter in subproblem_maxiter:
            # Solve Subproblem
            subprob = IterativeSubproblem(
                x=0,
                fun=lambda x: 0,
                jac=lambda x: J,
                hess=lambda x: H,
                k_easy=0.1,
                k_hard=0.2,
                maxiter=maxiter,
            )
            trust_radius = 1
            p, hits_boundary = subprob.solve(trust_radius)

            if maxiter is None:
                assert subprob.niter <= IterativeSubproblem.MAXITER_DEFAULT
            else:
                assert subprob.niter <= maxiter
