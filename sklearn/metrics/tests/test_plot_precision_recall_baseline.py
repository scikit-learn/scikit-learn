import matplotlib
import numpy as np

matplotlib.use("agg")
import matplotlib.pyplot as plt

from sklearn.metrics import PrecisionRecallDisplay
from sklearn.utils._testing import assert_allclose


def test_baseline_is_pos_ratio():
    y_true = np.array([1, 0, 0, 1, 0, 0, 0, 1, 0, 0])
    y_scores = np.linspace(0, 1, 10)

    disp = PrecisionRecallDisplay.from_predictions(y_true, y_scores)
    fig, ax = plt.subplots()
    disp.plot(ax=ax, plot_chance_level=True)

    lines = [line for line in ax.get_lines()]
    baseline_lines = [ln for ln in lines if ln.get_label().startswith("baseline")]
    assert len(baseline_lines) == 1, "Expected exactly one baseline line to be plotted"

    baseline_y = baseline_lines[0].get_ydata()[0]
    expected = float(np.count_nonzero(y_true == 1)) / float(len(y_true))
    assert_allclose(baseline_y, expected, rtol=1e-7, atol=1e-7)
