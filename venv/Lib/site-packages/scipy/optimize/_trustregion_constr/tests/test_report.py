from scipy.optimize import minimize, Bounds

def test_gh10880():
    # checks that verbose reporting works with trust-constr
    bnds = Bounds(1, 2)
    opts = {'maxiter': 1000, 'verbose': 2}
    minimize(lambda x: x**2, x0=2., method='trust-constr', bounds=bnds, options=opts)

    opts = {'maxiter': 1000, 'verbose': 3}
    minimize(lambda x: x**2, x0=2., method='trust-constr', bounds=bnds, options=opts)
