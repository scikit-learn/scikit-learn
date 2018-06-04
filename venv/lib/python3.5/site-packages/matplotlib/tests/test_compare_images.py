from __future__ import absolute_import, division, print_function

import six

import io
import os
import shutil
import warnings

from numpy.testing import assert_almost_equal
import pytest
from pytest import approx

from matplotlib.testing.compare import compare_images
from matplotlib.testing.decorators import _image_directories, image_comparison
from matplotlib.testing.exceptions import ImageComparisonFailure


baseline_dir, result_dir = _image_directories(lambda: 'dummy func')


# Tests of the image comparison algorithm.
@pytest.mark.parametrize(
    'im1, im2, tol, expect_rms',
    [
        # Comparison of an image and the same image with minor differences.
        # This expects the images to compare equal under normal tolerance, and
        # have a small RMS.
        ('basn3p02.png', 'basn3p02-minorchange.png', 10, None),
        # Now test with no tolerance.
        ('basn3p02.png', 'basn3p02-minorchange.png', 0, 6.50646),
        # Comparison with an image that is shifted by 1px in the X axis.
        ('basn3p02.png', 'basn3p02-1px-offset.png', 0, 90.15611),
        # Comparison with an image with half the pixels shifted by 1px in the X
        # axis.
        ('basn3p02.png', 'basn3p02-half-1px-offset.png', 0, 63.75),
        # Comparison of an image and the same image scrambled.
        # This expects the images to compare completely different, with a very
        # large RMS.
        # Note: The image has been scrambled in a specific way, by having
        # each color component of each pixel randomly placed somewhere in the
        # image. It contains exactly the same number of pixels of each color
        # value of R, G and B, but in a totally different position.
        # Test with no tolerance to make sure that we pick up even a very small
        # RMS error.
        ('basn3p02.png', 'basn3p02-scrambled.png', 0, 172.63582),
        # Comparison of an image and a slightly brighter image.
        # The two images are solid color, with the second image being exactly 1
        # color value brighter.
        # This expects the images to compare equal under normal tolerance, and
        # have an RMS of exactly 1.
        ('all127.png', 'all128.png', 0, 1),
        # Now test the reverse comparison.
        ('all128.png', 'all127.png', 0, 1),
    ])
def test_image_comparison_expect_rms(im1, im2, tol, expect_rms):
    """Compare two images, expecting a particular RMS error.

    im1 and im2 are filenames relative to the baseline_dir directory.

    tol is the tolerance to pass to compare_images.

    expect_rms is the expected RMS value, or None. If None, the test will
    succeed if compare_images succeeds. Otherwise, the test will succeed if
    compare_images fails and returns an RMS error almost equal to this value.
    """
    im1 = os.path.join(baseline_dir, im1)
    im2_src = os.path.join(baseline_dir, im2)
    im2 = os.path.join(result_dir, im2)
    # Move im2 from baseline_dir to result_dir. This will ensure that
    # compare_images writes the diff file to result_dir, instead of trying to
    # write to the (possibly read-only) baseline_dir.
    shutil.copyfile(im2_src, im2)
    results = compare_images(im1, im2, tol=tol, in_decorator=True)

    if expect_rms is None:
        assert results is None
    else:
        assert results is not None
        assert results['rms'] == approx(expect_rms, abs=1e-4)


# The following tests are used by test_nose_image_comparison to ensure that the
# image_comparison decorator continues to work with nose. They should not be
# prefixed by test_ so they don't run with pytest.


def nosetest_empty():
    pass


def nosetest_simple_figure():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4), dpi=100)
    ax.plot([1, 2, 3], [3, 4, 5])
    return fig


def nosetest_manual_text_removal():
    from matplotlib.testing.decorators import ImageComparisonTest

    fig = nosetest_simple_figure()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        # Make sure this removes text like it should.
        ImageComparisonTest.remove_text(fig)

    assert len(w) == 1
    assert 'remove_text function was deprecated in version 2.1.' in str(w[0])


@pytest.mark.parametrize(
    'func, kwargs, errors, failures, dots',
    [
        (nosetest_empty, {'baseline_images': []}, [], [], ''),
        (nosetest_empty, {'baseline_images': ['foo']},
         [(AssertionError,
           'Test generated 0 images but there are 1 baseline images')],
         [],
         'E'),
        (nosetest_simple_figure,
         {'baseline_images': ['basn3p02'], 'extensions': ['png'],
          'remove_text': True},
         [],
         [(ImageComparisonFailure, 'Image sizes do not match expected size:')],
         'F'),
        (nosetest_simple_figure,
         {'baseline_images': ['simple']},
         [],
         [(ImageComparisonFailure, 'images not close')] * 3,
         'FFF'),
        (nosetest_simple_figure,
         {'baseline_images': ['simple'], 'remove_text': True},
         [],
         [],
         '...'),
        (nosetest_manual_text_removal,
         {'baseline_images': ['simple']},
         [],
         [],
         '...'),
    ],
    ids=[
        'empty',
        'extra baselines',
        'incorrect shape',
        'failing figure',
        'passing figure',
        'manual text removal',
    ])
def test_nose_image_comparison(func, kwargs, errors, failures, dots,
                               monkeypatch):
    nose = pytest.importorskip('nose')
    monkeypatch.setattr('matplotlib._called_from_pytest', False)

    class TestResultVerifier(nose.result.TextTestResult):
        def __init__(self, *args, **kwargs):
            super(TestResultVerifier, self).__init__(*args, **kwargs)
            self.error_count = 0
            self.failure_count = 0

        def addError(self, test, err):
            super(TestResultVerifier, self).addError(test, err)

            if self.error_count < len(errors):
                assert err[0] is errors[self.error_count][0]
                assert errors[self.error_count][1] in str(err[1])
            else:
                raise err[1]
            self.error_count += 1

        def addFailure(self, test, err):
            super(TestResultVerifier, self).addFailure(test, err)

            assert self.failure_count < len(failures), err[1]
            assert err[0] is failures[self.failure_count][0]
            assert failures[self.failure_count][1] in str(err[1])
            self.failure_count += 1

    # Make sure that multiple extensions work, but don't require LaTeX or
    # Inkscape to do so.
    kwargs.setdefault('extensions', ['png', 'png', 'png'])

    func = image_comparison(**kwargs)(func)
    loader = nose.loader.TestLoader()
    suite = loader.loadTestsFromGenerator(
        func,
        'matplotlib.tests.test_compare_images')
    if six.PY2:
        output = io.BytesIO()
    else:
        output = io.StringIO()
    result = TestResultVerifier(stream=output, descriptions=True, verbosity=1)
    with warnings.catch_warnings():
        # Nose uses deprecated stuff; we don't care about it.
        warnings.simplefilter('ignore', DeprecationWarning)
        suite.run(result=result)

    assert output.getvalue() == dots
    assert result.error_count == len(errors)
    assert result.failure_count == len(failures)
