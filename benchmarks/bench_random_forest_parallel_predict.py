"""
This benchmark compares the predictions of a RandomForestRegressor
when the parallelization is disabled. The benchmark depends on
four parameters:

* number of estimators (-e)
* number of observations in the batch (-n)
* number of features (-f)
* number of repetitions (-r)
* use assume_finite (-a)
* onnx comparison (-o)

Option `-o` uses module *onnxruntime* and *mlprodict*
which provide two C++ implementation of the random forest
predict method. It requires a conversion of the model
into *onnx* done by modules *sklearn-onnx* and *onnx*.
The first execution fails after saving a trained model
to make sure the training is not part of the benchmark.
*py-psy* can be run using the following command line:

::

    py-spy record --native --function --rate=10 -o profile.svg --
    python bench_random_forest_parallel_predict.py -e 100 -n 1 -f 10 -r 1000
"""
import argparse
import time
import pickle
import os
import numpy as np
from sklearn import config_context
from sklearn.datasets import make_regression
from sklearn.ensemble import RandomForestRegressor


def build_model(e, n, f, o):
    filename = "RandomForestRegressor-e%d-n%d-f%d-onnx%d.pkl" % (e, n, f, o)
    if os.path.exists(filename):
        with open(filename, "rb") as f:
            return pickle.load(f)

    print("training e=%d n=%d f=%d onnx=%d" % (e, n, f, o))
    nt = 10000
    X, y = make_regression(
        n_samples=nt + n, n_features=f, n_informative=f // 2,
        n_targets=1, random_state=1)
    X_train, X_test = X[:nt], X[nt:]
    y_train = y[:nt]
    rf = RandomForestRegressor(n_estimators=e, random_state=1)
    rf.fit(X_train, y_train)

    data = dict(model=rf, data=X_test.astype(np.float32))
    if o:
        # compares with onnx
        print("convert to onnx")
        from skl2onnx import to_onnx
        model_onnx = to_onnx(rf, X_train[:1].astype(np.float32))
        buffer_onnx = model_onnx.SerializeToString()
        data['onnx'] = buffer_onnx

    print("saving to '%s'" % filename)
    with open(filename, "wb") as f:
        pickle.dump(data, f)
    raise RuntimeError("Data and data are cached. Run the script "
                       "again to measure the performance.")


def _run_predict(model, X, repeat, parallel_predict, check_input):
    for r in range(repeat):
        model.predict(X, parallel_predict=parallel_predict,
                      check_input=check_input)


def pp_true(model, X, repeat):
    _run_predict(model, X, repeat, True, True)


def pp_false(model, X, repeat):
    _run_predict(model, X, repeat, False, True)


def pp_false_nocheck(model, X, repeat):
    _run_predict(model, X, repeat, False, False)


def onnx_predict(sess, X, repeat):
    for r in range(repeat):
        sess.run(None, {'X': X})


def mlprodict_predict(sess, X, repeat):
    for r in range(repeat):
        sess.run({'X': X})


def benchmark(model, model_onnx, X, repeat):
    begin = time.perf_counter()
    pp_true(model, X, repeat)
    end1 = time.perf_counter()
    pp_false(model, X, repeat)
    end2 = time.perf_counter()
    pp_false_nocheck(model, X, repeat)
    end3 = time.perf_counter()
    r = repeat
    e = len(model.estimators_)
    n = X.shape[0]
    f = X.shape[1]
    print("parallel=T check=T:  r=%d e=%d n=%d f=%d time=%f" % (
        r, e, n, f, end1 - begin))
    print("parallel=F check=T:  r=%d e=%d n=%d f=%d time=%f" % (
        r, e, n, f, end2 - end1))
    print("parallel=F check=F:  r=%d e=%d n=%d f=%d time=%f" % (
        r, e, n, f, end3 - end2))

    if model_onnx is not None:
        from onnxruntime import InferenceSession
        from mlprodict.onnxrt import OnnxInference
        sess = InferenceSession(model_onnx)
        oinf = OnnxInference(model_onnx, runtime='python_compiled')

        begin = time.perf_counter()
        onnx_predict(sess, X, repeat)
        end = time.perf_counter()
        print("onnxruntime predict: r=%d e=%d n=%d f=%d time=%f" % (
            r, e, n, f, end - begin))

        begin = time.perf_counter()
        mlprodict_predict(oinf, X, repeat)
        end = time.perf_counter()
        print("mlprodict_predict:   r=%d e=%d n=%d f=%d time=%f" % (
            r, e, n, f, end - begin))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', type=int, default=200)
    parser.add_argument('-n', type=int, default=10000)
    parser.add_argument('-f', type=int, default=200)
    parser.add_argument('-r', type=int, default=10)
    parser.add_argument('-a', type=int, default=1)
    parser.add_argument('-o', type=int, default=1)
    args = parser.parse_args()

    model_data = build_model(args.e, args.n, args.f, args.o)

    if args.a:
        with config_context(assume_finite=True):
            benchmark(model_data['model'], model_data.get('onnx', None),
                      model_data['data'], args.r)
    else:
        benchmark(model_data['model'], model_data.get('onnx', None),
                  model_data['data'], args.r)


main()
