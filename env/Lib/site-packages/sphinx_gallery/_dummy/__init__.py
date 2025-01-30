from .nested import NestedDummyClass  # noqa: F401


class DummyClass:
    """Dummy class for testing method resolution."""

    def run(self):
        """Do nothing."""
        pass

    @property
    def prop(self):
        """Property."""
        return "Property"
