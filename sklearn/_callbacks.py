# License: BSD 3 clause
from typing import List, Callable
from abc import ABC, abstractmethod


def _eval_callbacks(callbacks: List[Callable], **kwargs) -> None:
    if callbacks is None:
        return

    for callback in callbacks:
        callback(**kwargs)


class BaseCallback(ABC):
    @abstractmethod
    def fit(self, X, y) -> None:
        pass

    @abstractmethod
    def __call__(self, **kwargs) -> None:
        pass
