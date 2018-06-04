from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import pytest

import matplotlib


def pytest_configure(config):
    matplotlib.use('agg')
    matplotlib._called_from_pytest = True
    matplotlib._init_tests()


def pytest_unconfigure(config):
    matplotlib._called_from_pytest = False


@pytest.fixture(autouse=True)
def mpl_test_settings(request):
    from matplotlib.testing.decorators import _do_cleanup

    original_units_registry = matplotlib.units.registry.copy()
    original_settings = matplotlib.rcParams.copy()

    backend = None
    backend_marker = request.keywords.get('backend')
    if backend_marker is not None:
        assert len(backend_marker.args) == 1, \
            "Marker 'backend' must specify 1 backend."
        backend = backend_marker.args[0]
        prev_backend = matplotlib.get_backend()

    style = '_classic_test'  # Default of cleanup and image_comparison too.
    style_marker = request.keywords.get('style')
    if style_marker is not None:
        assert len(style_marker.args) == 1, \
            "Marker 'style' must specify 1 style."
        style = style_marker.args[0]

    matplotlib.testing.setup()
    if backend is not None:
        # This import must come after setup() so it doesn't load the default
        # backend prematurely.
        import matplotlib.pyplot as plt
        plt.switch_backend(backend)
    matplotlib.style.use(style)
    try:
        yield
    finally:
        if backend is not None:
            plt.switch_backend(prev_backend)
        _do_cleanup(original_units_registry,
                    original_settings)


@pytest.fixture
def mpl_image_comparison_parameters(request, extension):
    # This fixture is applied automatically by the image_comparison decorator.
    #
    # The sole purpose of this fixture is to provide an indirect method of
    # obtaining parameters *without* modifying the decorated function
    # signature. In this way, the function signature can stay the same and
    # pytest won't get confused.
    # We annotate the decorated function with any parameters captured by this
    # fixture so that they can be used by the wrapper in image_comparison.
    baseline_images = request.keywords['baseline_images'].args[0]
    if baseline_images is None:
        # Allow baseline image list to be produced on the fly based on current
        # parametrization.
        baseline_images = request.getfixturevalue('baseline_images')

    func = request.function
    func.__wrapped__.parameters = (baseline_images, extension)
    try:
        yield
    finally:
        delattr(func.__wrapped__, 'parameters')


@pytest.fixture
def pd():
    """Fixture to import and configure pandas."""
    pd = pytest.importorskip('pandas')
    try:
        from pandas.plotting import (
            register_matplotlib_converters as register)
    except ImportError:
        from pandas.tseries.converter import register
    register()
    try:
        yield pd
    finally:
        try:
            from pandas.plotting import (
                deregister_matplotlib_converters as deregister)
        except ImportError:
            pass
        else:
            deregister()
