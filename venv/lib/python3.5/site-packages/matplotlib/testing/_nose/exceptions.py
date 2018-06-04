class KnownFailureTest(Exception):
    """
    Raise this exception to mark a test as a known failing test.
    """


class KnownFailureDidNotFailTest(Exception):
    """
    Raise this exception to mark a test should have failed but did not.
    """
