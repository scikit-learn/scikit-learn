import pytest
import warnings
import numpy as np


class Wrapper:
    def __init__(self, array):
        self.array = array

    def __len__(self):
        return len(self.array)

    def __getitem__(self, item):
        return type(self)(self.array[item])

    def __getattr__(self, name):
        if name.startswith("__array_"):
            warnings.warn("object got converted", UserWarning, stacklevel=1)

        return getattr(self.array, name)

    def __repr__(self):
        return "<Wrapper({self.array})>".format(self=self)

@pytest.mark.filterwarnings("error")
def test_getattr_warning():
    array = Wrapper(np.arange(10))
    with pytest.raises(UserWarning, match="object got converted"):
        np.asarray(array)
