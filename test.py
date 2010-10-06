# Small script to run the test suite.
#
# This is here for compatibility reasons: some tests depend on
# multiprocessing and thus must be run from a __main__ block,
# otherwise they will fail on windows systems

import nose

if __name__ == '__main__':
    try:
        import scikits.learn as skl
        skl.test()
    except ImportError:
        nose.run('scikits.learn')
