"""Test utilities for docstring."""

import pytest

from sklearn.utils._docstring import inject_docstring

expected_func_docstring = """A function.

        Parameters
        ----------
        xxx

        yyy"""

expected_cls_docstring = """A class.

        Parameters
        ----------
        xxx

        yyy"""


@pytest.fixture(scope="function")
def func():
    def func_to_inject(param_1, param_2):
        """A function.

        Parameters
        ----------
        {param_1}

        {param_2}
        """
        return param_1, param_2
    return func_to_inject, expected_func_docstring


@pytest.fixture(scope="function")
def Klass():
    class KlassToInject:
        """A class.

        Parameters
        ----------
        {param_1}

        {param_2}
        """
        def __init__(self, param_1, param_2):
            self.param_1 = param_1
            self.param_2 = param_2

    return KlassToInject, expected_cls_docstring


@pytest.fixture(
    params=["func", "Klass"],
    scope="function",
)
def object_to_inject(request):
    return request.getfixturevalue(request.param)


def test_inject_docstring(object_to_inject):
    obj, expected_docstring = object_to_inject
    obj_injected = inject_docstring(param_1="xxx", param_2="yyy")(
        obj
    )
    assert (obj_injected.__doc__.rstrip() ==
            expected_docstring.rstrip())
