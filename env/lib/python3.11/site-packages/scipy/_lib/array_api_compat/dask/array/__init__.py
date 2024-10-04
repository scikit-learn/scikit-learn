from dask.array import * # noqa: F403

# These imports may overwrite names from the import * above.
from ._aliases import * # noqa: F403

__array_api_version__ = '2022.12'

__import__(__package__ + '.linalg')
