import numpy as np

from sklearn.metrics._classification import batch_f1_score, f1_score


def test_batch_f1_score():
    """Test correctness of batch_f1_score."""
    y_true = np.array([0, 1, 1, 0, 1, 0])
    y_pred = np.array([0, 1, 0, 0, 1, 1])
    batch_score = batch_f1_score(y_true, y_pred, batch_size=2, average="macro")
    regular_score = f1_score(y_true, y_pred, average="macro")
    assert np.isclose(batch_score, regular_score, rtol=1e-05), (
        "Batch F1 score does not match regular F1 score"
    )
