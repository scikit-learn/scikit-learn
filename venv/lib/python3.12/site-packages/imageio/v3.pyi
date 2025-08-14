from typing import Any, Dict, Iterator, List, Literal, Optional, Union, overload

import numpy as np

from .core.imopen import imopen as imopen
from .core.v3_plugin_api import ImageProperties
from .typing import ArrayLike, ImageResource

def imread(
    uri: ImageResource,
    *,
    index: Optional[int] = 0,
    plugin: str = None,
    extension: str = None,
    format_hint: str = None,
    **kwargs,
) -> np.ndarray: ...
def imiter(
    uri: ImageResource,
    *,
    plugin: str = None,
    extension: str = None,
    format_hint: str = None,
    **kwargs,
) -> Iterator[np.ndarray]: ...
@overload
def imwrite(
    uri: Literal["<bytes>"],
    image: Union[ArrayLike, List[ArrayLike]],
    *,
    plugin: str = None,
    extension: str = None,
    format_hint: str = None,
    **kwargs,
) -> bytes: ...
@overload
def imwrite(
    uri: ImageResource,
    image: Union[ArrayLike, List[ArrayLike]],
    *,
    plugin: str = None,
    extension: str = None,
    format_hint: str = None,
    **kwargs,
) -> None: ...
def improps(
    uri,
    *,
    index: Optional[int] = 0,
    plugin: str = None,
    extension: str = None,
    **kwargs,
) -> ImageProperties: ...
def immeta(
    uri,
    *,
    index: Optional[int] = 0,
    plugin: str = None,
    extension: str = None,
    exclude_applied: bool = True,
    **kwargs,
) -> Dict[str, Any]: ...
