from __future__ import annotations

from collections.abc import Awaitable
from typing import TYPE_CHECKING, Any, Generic, TypeVar

from polars._utils.wrap import wrap_df
from polars.dependencies import _GEVENT_AVAILABLE

if TYPE_CHECKING:
    from asyncio.futures import Future
    from collections.abc import Generator

    from polars.polars import PyDataFrame


T = TypeVar("T")


class _GeventDataFrameResult(Generic[T]):
    __slots__ = ("_result", "_value", "_watcher")

    def __init__(self) -> None:
        if not _GEVENT_AVAILABLE:
            msg = (
                "gevent is required for using LazyFrame.collect_async(gevent=True) or"
                "polars.collect_all_async(gevent=True)"
            )
            raise ImportError(msg)

        from gevent.event import AsyncResult  # type: ignore[import-untyped]
        from gevent.hub import get_hub  # type: ignore[import-untyped]

        self._value: None | Exception | PyDataFrame | list[PyDataFrame] = None
        self._result = AsyncResult()

        self._watcher = get_hub().loop.async_()
        self._watcher.start(self._watcher_callback)

    def get(
        self,
        block: bool = True,  # noqa: FBT001
        timeout: float | int | None = None,
    ) -> T:
        return self.result.get(block=block, timeout=timeout)

    @property
    def result(self) -> Any:
        # required if we did not made any switches and just want results later
        # with block=False and possibly without timeout
        if self._value is not None and not self._result.ready():
            self._watcher_callback()
        return self._result

    def _watcher_callback(self) -> None:
        if isinstance(self._value, Exception):
            self._result.set_exception(self._value)
        else:
            self._result.set(self._value)
        self._watcher.close()

    def _callback(self, obj: PyDataFrame | Exception) -> None:
        if not isinstance(obj, Exception):
            obj = wrap_df(obj)
        self._value = obj
        self._watcher.send()

    def _callback_all(self, obj: list[PyDataFrame] | Exception) -> None:
        if not isinstance(obj, Exception):
            obj = [wrap_df(pydf) for pydf in obj]
        self._value = obj
        self._watcher.send()


class _AioDataFrameResult(Awaitable[T], Generic[T]):
    __slots__ = ("loop", "result")

    def __init__(self) -> None:
        from asyncio import get_event_loop

        self.loop = get_event_loop()
        self.result: Future[T] = self.loop.create_future()

    def __await__(self) -> Generator[Any, None, T]:
        return self.result.__await__()

    def _callback(self, obj: PyDataFrame | Exception) -> None:
        if isinstance(obj, Exception):
            self.loop.call_soon_threadsafe(self.result.set_exception, obj)
        else:
            self.loop.call_soon_threadsafe(
                self.result.set_result,  # type: ignore[arg-type]
                wrap_df(obj),
            )

    def _callback_all(self, obj: list[PyDataFrame] | Exception) -> None:
        if isinstance(obj, Exception):
            self.loop.call_soon_threadsafe(self.result.set_exception, obj)
        else:
            self.loop.call_soon_threadsafe(
                self.result.set_result,  # type: ignore[arg-type]
                [wrap_df(pydf) for pydf in obj],
            )
