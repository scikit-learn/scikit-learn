from .._docstrings import DocstringComponents


EXAMPLE_DICT = dict(
    param_a="""
a : str
    The first parameter.
    """,
)


class ExampleClass:
    def example_method(self):
        """An example method.

        Parameters
        ----------
        a : str
           A method parameter.

        """


def example_func():
    """An example function.

    Parameters
    ----------
    a : str
        A function parameter.

    """


class TestDocstringComponents:

    def test_from_dict(self):

        obj = DocstringComponents(EXAMPLE_DICT)
        assert obj.param_a == "a : str\n    The first parameter."

    def test_from_nested_components(self):

        obj_inner = DocstringComponents(EXAMPLE_DICT)
        obj_outer = DocstringComponents.from_nested_components(inner=obj_inner)
        assert obj_outer.inner.param_a == "a : str\n    The first parameter."

    def test_from_function(self):

        obj = DocstringComponents.from_function_params(example_func)
        assert obj.a == "a : str\n    A function parameter."

    def test_from_method(self):

        obj = DocstringComponents.from_function_params(
            ExampleClass.example_method
        )
        assert obj.a == "a : str\n    A method parameter."
