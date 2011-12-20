"""
Helper functions for benchmarking
"""


def total_seconds(delta):
    """
    helper function to emulate function total_seconds,
    introduced in python2.7

    http://docs.python.org/library/datetime.html\
#datetime.timedelta.total_seconds
    """

    mu_sec = 1e-6  # number of seconds in one microseconds

    return delta.seconds + delta.microseconds * mu_sec
