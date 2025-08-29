import pytest

# List of files that pytest should ignore
collect_ignore = [
    "io/_plugins",
]


@pytest.fixture(autouse=True)
def handle_np2():
    # TODO: remove when we require numpy >= 2
    try:
        import numpy as np

        np.set_printoptions(legacy="1.21")
    except ImportError:
        pass
