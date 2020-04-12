# License: BSD 3 clause
from abc import ABC, abstractmethod


def _eval_callbacks(callbacks, **kwargs):
    if callbacks is None:
        return

    for callback in callbacks:
        callback(**kwargs)


class Callback(ABC):
    @abstractmethod
    def fit(self, X, y):
        pass

    @abstractmethod
    def __call__(self, **kwargs):
        pass


class ProgressBar(Callback):
    def __init__(self):
        self.pbar = None

    def fit(X, y):
        pass

    def __call__(self, **kwargs):
        from tqdm.auto import tqdm

        if self.pbar is None:
            self.pbar = tqdm(total=kwargs.get("n_iter_total"))
        self.pbar.update(1)
