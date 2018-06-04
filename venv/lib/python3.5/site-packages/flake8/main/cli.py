"""Command-line implementation of flake8."""
from flake8.main import application


def main(argv=None):
    # type: (Union[NoneType, List[str]]) -> NoneType
    """Execute the main bit of the application.

    This handles the creation of an instance of :class:`Application`, runs it,
    and then exits the application.

    :param list argv:
        The arguments to be passed to the application for parsing.
    """
    app = application.Application()
    app.run(argv)
    app.exit()
