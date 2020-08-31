from typing import TYPE_CHECKING
from typing import Union

import numpy as np

RandomState = Union[int, np.random.RandomState, None]

if TYPE_CHECKING:
    from typing_extensions import Literal  # type: ignore  # noqa
else:
    Literal = None
