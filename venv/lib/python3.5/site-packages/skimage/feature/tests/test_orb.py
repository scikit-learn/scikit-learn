import numpy as np
from skimage._shared.testing import assert_equal, assert_almost_equal
from skimage.feature import ORB
from skimage._shared import testing
from skimage import data
from skimage._shared.testing import test_parallel


img = data.coins()


@test_parallel()
def test_keypoints_orb_desired_no_of_keypoints():
    detector_extractor = ORB(n_keypoints=10, fast_n=12, fast_threshold=0.20)
    detector_extractor.detect(img)

    exp_rows = np.array([ 141.   ,  108.   ,  214.56 ,  131.   ,  214.272,
                           67.   ,  206.   ,  177.   ,  108.   ,  141.   ])
    exp_cols = np.array([ 323.   ,  328.   ,  282.24 ,  292.   ,  281.664,
                           85.   ,  260.   ,  284.   ,  328.8  ,  267.   ])

    exp_scales = np.array([1,  1,  1.44,  1,  1.728, 1, 1, 1, 1.2, 1])

    exp_orientations = np.array([ -53.97446153,   59.5055285 ,  -96.01885186,
                                 -149.70789506,  -94.70171899,  -45.76429535,
                                  -51.49752849,  113.57081195,   63.30428063,
                                  -79.56091118])
    exp_response = np.array([ 1.01168357,  0.82934145,  0.67784179,  0.57176438,
                              0.56637459,  0.52248355,  0.43696175,  0.42992376,
                              0.37700486,  0.36126832])

    assert_almost_equal(exp_rows, detector_extractor.keypoints[:, 0])
    assert_almost_equal(exp_cols, detector_extractor.keypoints[:, 1])
    assert_almost_equal(exp_scales, detector_extractor.scales)
    assert_almost_equal(exp_response, detector_extractor.responses)
    assert_almost_equal(exp_orientations,
                        np.rad2deg(detector_extractor.orientations), 5)

    detector_extractor.detect_and_extract(img)
    assert_almost_equal(exp_rows, detector_extractor.keypoints[:, 0])
    assert_almost_equal(exp_cols, detector_extractor.keypoints[:, 1])


def test_keypoints_orb_less_than_desired_no_of_keypoints():
    detector_extractor = ORB(n_keypoints=15, fast_n=12,
                             fast_threshold=0.33, downscale=2, n_scales=2)
    detector_extractor.detect(img)

    exp_rows = np.array([  58.,   65.,  108.,  140.,  203.])
    exp_cols = np.array([ 291.,  130.,  293.,  202.,  267.])

    exp_scales = np.array([1., 1., 1., 1., 1.])

    exp_orientations = np.array([-158.26941428,  -59.42996346,  151.93905955,
                                  -79.46341354,  -56.90052451])

    exp_response = np.array([ 0.2667641 ,  0.04009017, -0.17641695, -0.03243431,
                              0.26521259])

    assert_almost_equal(exp_rows, detector_extractor.keypoints[:, 0])
    assert_almost_equal(exp_cols, detector_extractor.keypoints[:, 1])
    assert_almost_equal(exp_scales, detector_extractor.scales)
    assert_almost_equal(exp_response, detector_extractor.responses)
    assert_almost_equal(exp_orientations,
                        np.rad2deg(detector_extractor.orientations), 5)

    detector_extractor.detect_and_extract(img)
    assert_almost_equal(exp_rows, detector_extractor.keypoints[:, 0])
    assert_almost_equal(exp_cols, detector_extractor.keypoints[:, 1])


def test_descriptor_orb():
    detector_extractor = ORB(fast_n=12, fast_threshold=0.20)

    exp_descriptors = np.array([[0, 1, 1, 1, 0, 1, 0, 1, 0, 1],
                                [1, 1, 1, 0, 0, 1, 0, 0, 1, 1],
                                [1, 0, 1, 1, 0, 0, 1, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0, 1, 0],
                                [1, 1, 0, 1, 1, 1, 0, 0, 1, 1],
                                [1, 1, 0, 1, 0, 0, 1, 0, 1, 1],
                                [0, 0, 1, 0, 1, 0, 0, 1, 1, 0],
                                [1, 0, 0, 0, 1, 0, 0, 0, 0, 1],
                                [0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                                [1, 1, 0, 1, 0, 1, 0, 0, 1, 1],
                                [1, 1, 1, 0, 0, 0, 1, 1, 1, 0],
                                [1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
                                [1, 1, 1, 0, 1, 1, 1, 1, 0, 0],
                                [1, 1, 0, 0, 1, 0, 0, 1, 0, 1],
                                [1, 1, 0, 0, 0, 0, 1, 0, 0, 1],
                                [0, 0, 0, 0, 1, 1, 1, 0, 1, 0],
                                [0, 0, 0, 0, 1, 1, 1, 0, 0, 1],
                                [0, 0, 0, 0, 0, 1, 1, 0, 1, 1],
                                [0, 0, 0, 0, 1, 0, 1, 0, 1, 1]], dtype=bool)
    detector_extractor.detect(img)
    detector_extractor.extract(img, detector_extractor.keypoints,
                               detector_extractor.scales,
                               detector_extractor.orientations)
    assert_equal(exp_descriptors,
                 detector_extractor.descriptors[100:120, 10:20])

    detector_extractor.detect_and_extract(img)
    assert_equal(exp_descriptors,
                 detector_extractor.descriptors[100:120, 10:20])


def test_no_descriptors_extracted_orb():
    img = np.ones((128, 128))
    detector_extractor = ORB()
    with testing.raises(RuntimeError):
        detector_extractor.detect_and_extract(img)
