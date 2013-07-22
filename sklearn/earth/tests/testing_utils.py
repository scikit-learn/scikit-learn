from functools import wraps
from nose import SkipTest

def if_statsmodels(func):
    """Test decorator that skips test if statsmodels not installed. """

    @wraps(func)
    def run_test(*args, **kwargs):
        try:
            import statsmodels
        except ImportError:
            raise SkipTest('statsmodels not available.')
        else:
            return func(*args, **kwargs)
    return run_test

def if_pandas(func):
    """Test decorator that skips test if pandas not installed. """

    @wraps(func)
    def run_test(*args, **kwargs):
        try:
            import pandas
        except ImportError:
            raise SkipTest('pandas not available.')
        else:
            return func(*args, **kwargs)
    return run_test

def if_patsy(func):
    """Test decorator that skips test if patsy not installed. """

    @wraps(func)
    def run_test(*args, **kwargs):
        try:
            import patsy
        except ImportError:
            raise SkipTest('patsy not available.')
        else:
            return func(*args, **kwargs)
    return run_test
