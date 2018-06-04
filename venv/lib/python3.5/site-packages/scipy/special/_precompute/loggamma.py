"""Precompute series coefficients for log-Gamma."""
from __future__ import division, print_function, absolute_import

try:
    import mpmath
except ImportError:
    pass


def stirling_series(N):
    coeffs = []
    with mpmath.workdps(100):
        for n in range(1, N + 1):
            coeffs.append(mpmath.bernoulli(2*n)/(2*n*(2*n - 1)))
    return coeffs


def taylor_series_at_1(N):
    coeffs = []
    with mpmath.workdps(100):
        coeffs.append(-mpmath.euler)
        for n in range(2, N + 1):
            coeffs.append((-1)**n*mpmath.zeta(n)/n)
    return coeffs


def main():
    print(__doc__)
    print()
    stirling_coeffs = [mpmath.nstr(x, 20, min_fixed=0, max_fixed=0)
                       for x in stirling_series(8)[::-1]]
    taylor_coeffs = [mpmath.nstr(x, 20, min_fixed=0, max_fixed=0)
                     for x in taylor_series_at_1(23)[::-1]]
    print("Stirling series coefficients")
    print("----------------------------")
    print("\n".join(stirling_coeffs))
    print()
    print("Taylor series coefficients")
    print("--------------------------")
    print("\n".join(taylor_coeffs))
    print()


if __name__ == '__main__':
    main()
    
