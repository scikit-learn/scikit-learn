import warnings

from pandas.util._exceptions import find_stack_level

from pandas._testing import *  # noqa:F401,F403,PDF014

warnings.warn(
    (
        "pandas.util.testing is deprecated. Use the functions in the "
        "public API at pandas.testing instead."
    ),
    FutureWarning,
    stacklevel=find_stack_level(),
)
