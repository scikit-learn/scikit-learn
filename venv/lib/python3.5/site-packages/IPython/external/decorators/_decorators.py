"""
Decorators for labeling and modifying behavior of test objects.

Decorators that merely return a modified version of the original
function object are straightforward. Decorators that return a new
function object need to use
::

  nose.tools.make_decorator(original_function)(decorator)

in returning the decorator, in order to preserve meta-data such as
function name, setup and teardown functions and so on - see
``nose.tools`` for more information.

"""

# IPython changes: make this work if numpy not available
# Original code:
try:
    from ._numpy_testing_noseclasses import KnownFailureTest
except:
    pass

# End IPython changes


def skipif(skip_condition, msg=None):
    """
    Make function raise SkipTest exception if a given condition is true.

    If the condition is a callable, it is used at runtime to dynamically
    make the decision. This is useful for tests that may require costly
    imports, to delay the cost until the test suite is actually executed.

    Parameters
    ----------
    skip_condition : bool or callable
        Flag to determine whether to skip the decorated test.
    msg : str, optional
        Message to give on raising a SkipTest exception. Default is None.

    Returns
    -------
    decorator : function
        Decorator which, when applied to a function, causes SkipTest
        to be raised when `skip_condition` is True, and the function
        to be called normally otherwise.

    Notes
    -----
    The decorator itself is decorated with the ``nose.tools.make_decorator``
    function in order to transmit function name, and various other metadata.

    """

    def skip_decorator(f):
        # Local import to avoid a hard nose dependency and only incur the
        # import time overhead at actual test-time.
        import nose

        # Allow for both boolean or callable skip conditions.
        if callable(skip_condition):
            skip_val = lambda : skip_condition()
        else:
            skip_val = lambda : skip_condition

        def get_msg(func,msg=None):
            """Skip message with information about function being skipped."""
            if msg is None:
                out = 'Test skipped due to test condition'
            else:
                out = '\n'+msg

            return "Skipping test: %s%s" % (func.__name__,out)

        # We need to define *two* skippers because Python doesn't allow both
        # return with value and yield inside the same function.
        def skipper_func(*args, **kwargs):
            """Skipper for normal test functions."""
            if skip_val():
                raise nose.SkipTest(get_msg(f,msg))
            else:
                return f(*args, **kwargs)

        def skipper_gen(*args, **kwargs):
            """Skipper for test generators."""
            if skip_val():
                raise nose.SkipTest(get_msg(f,msg))
            else:
                for x in f(*args, **kwargs):
                    yield x

        # Choose the right skipper to use when building the actual decorator.
        if nose.util.isgenerator(f):
            skipper = skipper_gen
        else:
            skipper = skipper_func

        return nose.tools.make_decorator(f)(skipper)

    return skip_decorator

def knownfailureif(fail_condition, msg=None):
    """
    Make function raise KnownFailureTest exception if given condition is true.

    Parameters
    ----------
    fail_condition : bool
        Flag to determine whether to mark the decorated test as a known
        failure (if True) or not (if False).
    msg : str, optional
        Message to give on raising a KnownFailureTest exception.
        Default is None.

    Returns
    -------
    decorator : function
        Decorator, which, when applied to a function, causes KnownFailureTest to
        be raised when `fail_condition` is True and the test fails.

    Notes
    -----
    The decorator itself is decorated with the ``nose.tools.make_decorator``
    function in order to transmit function name, and various other metadata.

    """
    if msg is None:
        msg = 'Test skipped due to known failure'

    def knownfail_decorator(f):
        # Local import to avoid a hard nose dependency and only incur the
        # import time overhead at actual test-time.
        import nose

        def knownfailer(*args, **kwargs):
            if fail_condition:
                raise KnownFailureTest(msg)
            else:
                return f(*args, **kwargs)
        return nose.tools.make_decorator(f)(knownfailer)

    return knownfail_decorator
