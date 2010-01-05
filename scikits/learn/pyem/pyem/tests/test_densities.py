def generate_test_data(n, d, mode = 'diag', file='test.dat'):
    """Generate a set of data of dimension d, with n frames,
    that is input data, mean, var and output of gden, so that
    other implementations can be tested against"""
    mu  = randn(1, d)
    if mode == 'diag':
        va  = abs(randn(1, d))
    elif mode == 'full':
        va  = randn(d, d)
        va  = dot(va, va.transpose())

    input   = randn(n, d)
    output  = gauss_den(input, mu, va)

    import tables
    h5file  = tables.openFile(file, "w")

    h5file.createArray(h5file.root, 'input', input)
    h5file.createArray(h5file.root, 'mu', mu)
    h5file.createArray(h5file.root, 'va', va)
    h5file.createArray(h5file.root, 'output', output)

    h5file.close()

def test_gauss_den():
    """"""
    # import tables
    # import numpy as N
    # 
    # filename    = 'dendata.h5'

    # # # Dimension 1
    # # d   = 1
    # # mu  = 1.0
    # # va  = 2.0

    # # X   = randn(1e3, 1)

    # # Y   = gauss_den(X, mu, va)

    # # h5file      = tables.openFile(filename, "w")

    # # h5file.createArray(h5file.root, 'X', X)
    # # h5file.createArray(h5file.root, 'mu', mu)
    # # h5file.createArray(h5file.root, 'va', va)
    # # h5file.createArray(h5file.root, 'Y', Y)

    # # h5file.close()

    # # # Dimension 2, diag
    # # d   = 2
    # # mu  = N.array([1.0, -2.0])
    # # va  = N.array([1.0, 2.0])

    # # X   = randn(1e3, 2)

    # # Y   = gauss_den(X, mu, va)

    # # h5file      = tables.openFile(filename, "w")

    # # h5file.createArray(h5file.root, 'X', X)
    # # h5file.createArray(h5file.root, 'mu', mu)
    # # h5file.createArray(h5file.root, 'va', va)
    # # h5file.createArray(h5file.root, 'Y', Y)

    # # Dimension 2, full
    # d   = 2
    # mu  = N.array([[0.2, -1.0]])
    # va  = N.array([[1.2, 0.1], [0.1, 0.5]])

    # X   = randn(1e3, 2)

    # Y   = gauss_den(X, mu, va)

    # h5file      = tables.openFile(filename, "w")

    # h5file.createArray(h5file.root, 'X', X)
    # h5file.createArray(h5file.root, 'mu', mu)
    # h5file.createArray(h5file.root, 'va', va)
    # h5file.createArray(h5file.root, 'Y', Y)

    # h5file.close()

    import numpy.testing as testing
    #=================
    # Small test in 1d
    #=================
    va  = 2.0
    mu  = 1.0
    X   = N.linspace(-2, 2, 10)[:, N.NewAxis]

    Yt  = N.array([0.02973257230591, 0.05512079811082, 0.09257745306945, 
            0.14086453882683,
            0.19418015562214, 0.24250166773127, 0.27436665745048, 0.28122547107069,
            0.26114678964743, 0.21969564473386])

    Y   = gauss_den(X, mu, va)
    try:
        testing.assert_array_almost_equal(Y, Yt)
        print "1d test succeded"
    except AssertionError:
        print "test fails in 1d"

    #============================
    # Small test in 2d (diagonal)
    #============================
    mu  = N.atleast_2d([-1.0, 2.0])
    va  = N.atleast_2d([2.0, 3.0])
    X1  = N.linspace(-2, 2, 10)[:, N.NewAxis]
    X2  = N.linspace(-1, 3, 10)[:, N.NewAxis]
    X   = N.concatenate(([X1, X2]), 1)
    
    Yt  = N.array([0.01129091565384, 0.02025416837152, 0.03081845516786, 
            0.03977576221540, 0.04354490552910, 0.04043592581117, 
            0.03184994053539, 0.02127948225225, 0.01205937178755, 
            0.00579694938623 ])

    Y   = gauss_den(X, mu, va)
    try:
        testing.assert_array_almost_equal(Y, Yt)
        print "2d diag test succeded"
    except AssertionError:
        print "test fails in 2d diag"

    #============================
    # Small test in 2d (full mat)
    #============================
    mu  = N.array([[0.2, -1.0]])
    va  = N.array([[1.2, 0.1], [0.1, 0.5]])
    X1  = N.linspace(-2, 2, 10)[:, N.NewAxis]
    X2  = N.linspace(-3, 3, 10)[:, N.NewAxis]
    X   = N.concatenate(([X1, X2]), 1)
    
    Yt  = N.array([0.00096157109751, 0.01368908714856,
        0.07380823191162, 0.15072050533842, 
        0.11656739937861, 0.03414436965525,
        0.00378789836599, 0.00015915297541, 
        0.00000253261067, 0.00000001526368])

    Y   = gauss_den(X, mu, va)
    try:
        testing.assert_array_almost_equal(Y, Yt)
        print "2d full test succeded"
    except AssertionError:
        print "test fails in 2d full"
