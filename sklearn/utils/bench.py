"""
Helper functions for benchmarking
"""


def total_seconds(delta):
    """
    helper function to emulate function total_seconds,
    introduced in python2.7

    https://docs.python.org/library/datetime.html\
#datetime.timedelta.total_seconds

    Parameters
    ----------
    delta : datetime object

    Returns
    -------
    int
        The number of seconds contained in delta
    """

    mu_sec = 1e-6  # number of seconds in one microseconds

    return delta.seconds + delta.microseconds * mu_sec
