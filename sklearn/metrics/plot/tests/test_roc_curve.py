import pytest


from sklearn.metrics import plot_roc_curve



def test_plot_roc_curve_error():
    bad_response_method = 'hello_world'

    plot_roc_curve()



def test_plot_roc_curve(pyplot):
    pass
