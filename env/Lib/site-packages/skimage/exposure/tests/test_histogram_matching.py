import numpy as np
import pytest
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from skimage import data
from skimage import exposure
from skimage._shared.utils import _supported_float_type
from skimage.exposure import histogram_matching


@pytest.mark.parametrize(
    'array, template, expected_array',
    [
        (np.arange(10), np.arange(100), np.arange(9, 100, 10)),
        (np.random.rand(4), np.ones(3), np.ones(4)),
    ],
)
def test_match_array_values(array, template, expected_array):
    # when
    matched = histogram_matching._match_cumulative_cdf(array, template)

    # then
    assert_array_almost_equal(matched, expected_array)


class TestMatchHistogram:
    image_rgb = data.chelsea()
    template_rgb = data.astronaut()

    @pytest.mark.parametrize(
        'image, reference, channel_axis',
        [
            (image_rgb, template_rgb, -1),
            (image_rgb[:, :, 0], template_rgb[:, :, 0], None),
        ],
    )
    def test_match_histograms(self, image, reference, channel_axis):
        """Assert that pdf of matched image is close to the reference's pdf for
        all channels and all values of matched"""

        matched = exposure.match_histograms(image, reference, channel_axis=channel_axis)

        matched_pdf = self._calculate_image_empirical_pdf(matched)
        reference_pdf = self._calculate_image_empirical_pdf(reference)

        for channel in range(len(matched_pdf)):
            reference_values, reference_quantiles = reference_pdf[channel]
            matched_values, matched_quantiles = matched_pdf[channel]

            for i, matched_value in enumerate(matched_values):
                closest_id = (np.abs(reference_values - matched_value)).argmin()
                assert_almost_equal(
                    matched_quantiles[i], reference_quantiles[closest_id], decimal=1
                )

    @pytest.mark.parametrize('channel_axis', (0, 1, -1))
    def test_match_histograms_channel_axis(self, channel_axis):
        """Assert that pdf of matched image is close to the reference's pdf for
        all channels and all values of matched"""

        image = np.moveaxis(self.image_rgb, -1, channel_axis)
        reference = np.moveaxis(self.template_rgb, -1, channel_axis)
        matched = exposure.match_histograms(image, reference, channel_axis=channel_axis)
        assert matched.dtype == image.dtype
        matched = np.moveaxis(matched, channel_axis, -1)
        reference = np.moveaxis(reference, channel_axis, -1)
        matched_pdf = self._calculate_image_empirical_pdf(matched)
        reference_pdf = self._calculate_image_empirical_pdf(reference)

        for channel in range(len(matched_pdf)):
            reference_values, reference_quantiles = reference_pdf[channel]
            matched_values, matched_quantiles = matched_pdf[channel]

            for i, matched_value in enumerate(matched_values):
                closest_id = (np.abs(reference_values - matched_value)).argmin()
                assert_almost_equal(
                    matched_quantiles[i], reference_quantiles[closest_id], decimal=1
                )

    @pytest.mark.parametrize('dtype', [np.float16, np.float32, np.float64])
    def test_match_histograms_float_dtype(self, dtype):
        """float16 or float32 inputs give float32 output"""
        image = self.image_rgb.astype(dtype, copy=False)
        reference = self.template_rgb.astype(dtype, copy=False)
        matched = exposure.match_histograms(image, reference)
        assert matched.dtype == _supported_float_type(dtype)

    @pytest.mark.parametrize(
        'image, reference',
        [(image_rgb, template_rgb[:, :, 0]), (image_rgb[:, :, 0], template_rgb)],
    )
    def test_raises_value_error_on_channels_mismatch(self, image, reference):
        with pytest.raises(ValueError):
            exposure.match_histograms(image, reference)

    @classmethod
    def _calculate_image_empirical_pdf(cls, image):
        """Helper function for calculating empirical probability density
        function of a given image for all channels"""

        if image.ndim > 2:
            image = image.transpose(2, 0, 1)
        channels = np.array(image, copy=False, ndmin=3)

        channels_pdf = []
        for channel in channels:
            channel_values, counts = np.unique(channel, return_counts=True)
            channel_quantiles = np.cumsum(counts).astype(np.float64)
            channel_quantiles /= channel_quantiles[-1]

            channels_pdf.append((channel_values, channel_quantiles))

        return np.asarray(channels_pdf, dtype=object)

    def test_match_histograms_consistency(self):
        """ensure equivalent results for float and integer-based code paths"""
        image_u8 = self.image_rgb
        reference_u8 = self.template_rgb
        image_f64 = self.image_rgb.astype(np.float64)
        reference_f64 = self.template_rgb.astype(np.float64, copy=False)
        matched_u8 = exposure.match_histograms(image_u8, reference_u8)
        matched_f64 = exposure.match_histograms(image_f64, reference_f64)
        assert_array_almost_equal(matched_u8.astype(np.float64), matched_f64)
