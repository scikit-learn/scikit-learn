from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Mapping

    from rmm.mr import DeviceMemoryResource  # type: ignore[import-not-found]


class GPUEngine:
    """
    Configuration options for the GPU execution engine.

    Use this if you want control over details of the execution.

    Parameters
    ----------
    device : int, default None
        Select the GPU used to run the query. If not provided, the
        query uses the current CUDA device.
    memory_resource : rmm.mr.DeviceMemoryResource, default None
        Provide a memory resource for GPU memory allocations.

        .. warning::
           If passing a `memory_resource`, you must ensure that it is valid
           for the selected `device`. See the `RMM documentation
           <https://github.com/rapidsai/rmm?tab=readme-ov-file#multiple-devices>`_
           for more details.

    raise_on_fail : bool, default False
        If True, do not fall back to the Polars CPU engine if the GPU
        engine cannot execute the query, but instead raise an error.

    """

    device: int | None
    """Device on which to run query."""
    memory_resource: DeviceMemoryResource | None
    """Memory resource to use for device allocations."""
    raise_on_fail: bool
    """
    Whether unsupported queries should raise an error, rather than falling
    back to the CPU engine.
    """
    config: Mapping[str, Any]
    """Additional configuration options for the engine."""

    def __init__(
        self,
        *,
        device: int | None = None,
        memory_resource: Any | None = None,
        raise_on_fail: bool = False,
        **kwargs: Any,
    ) -> None:
        self.device = device
        self.memory_resource = memory_resource
        # Avoids need for changes in cudf-polars
        kwargs["raise_on_fail"] = raise_on_fail
        self.config = kwargs
