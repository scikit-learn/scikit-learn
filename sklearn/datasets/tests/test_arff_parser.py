from sklearn.datasets import fetch_openml


def test_xxx():
    # df = fetch_openml("iris", version=1, as_frame=True, parser="pandas")
    # df = fetch_openml("mnist_784", version=1, as_frame=False, parser="liac-arff")
    df = fetch_openml("mnist_784", version=1, as_frame=True, parser="pandas")
