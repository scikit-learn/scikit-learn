from sklearn.datasets import fetch_openml


def test_xxx():
    # df = fetch_openml("titanic", version=1, as_frame=False, parser="liac-arff")
    # df = fetch_openml("mnist_784", version=1, as_frame=False, parser="liac-arff")
    df = fetch_openml("mnist_784", version=1, as_frame=False, parser="pandas")
