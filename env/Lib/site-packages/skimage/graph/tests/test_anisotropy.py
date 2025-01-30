import numpy as np
import skimage.graph.mcp as mcp

from skimage._shared.testing import assert_array_equal


a = np.ones((8, 8), dtype=np.float32)

horizontal_ramp = np.array(
    [
        [
            0.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
        ],
        [
            0.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
        ],
        [
            0.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
        ],
        [
            0.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
        ],
        [
            0.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
        ],
        [
            0.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
        ],
        [
            0.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
        ],
        [
            0.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
        ],
    ]
)

vertical_ramp = np.array(
    [
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        [
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
        ],
        [
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
        ],
        [
            3.0,
            3.0,
            3.0,
            3.0,
            3.0,
            3.0,
            3.0,
            3.0,
        ],
        [
            4.0,
            4.0,
            4.0,
            4.0,
            4.0,
            4.0,
            4.0,
            4.0,
        ],
        [
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
        ],
        [
            6.0,
            6.0,
            6.0,
            6.0,
            6.0,
            6.0,
            6.0,
            6.0,
        ],
        [
            7.0,
            7.0,
            7.0,
            7.0,
            7.0,
            7.0,
            7.0,
            7.0,
        ],
    ]
)


def test_anisotropy():
    # Create seeds; vertical seeds create a horizonral ramp
    seeds_for_horizontal = [(i, 0) for i in range(8)]
    seeds_for_vertcal = [(0, i) for i in range(8)]

    for sy in range(1, 5):
        for sx in range(1, 5):
            sampling = sy, sx
            # Trace horizontally
            m1 = mcp.MCP_Geometric(a, sampling=sampling, fully_connected=True)
            costs1, traceback = m1.find_costs(seeds_for_horizontal)
            # Trace vertically
            m2 = mcp.MCP_Geometric(a, sampling=sampling, fully_connected=True)
            costs2, traceback = m2.find_costs(seeds_for_vertcal)

            # Check
            assert_array_equal(costs1, horizontal_ramp * sx)
            assert_array_equal(costs2, vertical_ramp * sy)
