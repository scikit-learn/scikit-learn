import pytest
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay

def test_plot_sets_title():
    fpr = np.array([0.0, 0.1, 0.2, 1.0])
    tpr = np.array([0.0, 0.4, 0.8, 1.0])
    display = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=0.85, name="Test")

    fig, ax = plt.subplots()
    title_text = "Custom ROC Title"
    display.plot(ax=ax, title=title_text)

    assert ax.get_title() == title_text
