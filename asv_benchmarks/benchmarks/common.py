import os
import json
import timeit
import pickle
import itertools
from multiprocessing import cpu_count
from abc import ABC, abstractmethod

import numpy as np


def get_from_config():
    current_path = os.path.dirname(os.path.realpath(__file__))

    config_path = os.path.join(current_path, 'config.json')
    with open(config_path, 'r') as config_file:
        config_file = ''.join(line for line in config_file
                              if line and '//' not in line)
        config = json.loads(config_file)

    profile = config['profile']

    n_jobs_vals = config['n_jobs_vals']
    if not n_jobs_vals:
        n_jobs_vals = list(range(1, 1 + cpu_count()))

    cache_path = os.path.join(current_path, 'cache')
    if not os.path.exists(cache_path):
        os.mkdir(cache_path)
    estimators_path = os.path.join(current_path, 'cache', 'estimators')
    if not os.path.exists(estimators_path):
        os.mkdir(estimators_path)
    tmp_path = os.path.join(current_path, 'cache', 'tmp')
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)

    save_estimators = config['save_estimators']
    save_folder = os.getenv('ASV_COMMIT', 'new')[:8]

    if save_estimators:
        save_path = os.path.join(estimators_path, save_folder)
        if not os.path.exists(save_path):
            os.mkdir(save_path)

    base_folder = config['base_folder']

    bench_predict = config['bench_predict']
    bench_transform = config['bench_transform']

    return (profile, n_jobs_vals, save_estimators, save_folder, base_folder,
            bench_predict, bench_transform)


def get_estimator_path(benchmark, folder, params, save=False):
    folder = os.path.join('estimators', folder) if save else 'tmp'
    f_name = (benchmark.__class__.__name__[:-6]
              + '_estimator_' + '_'.join(list(map(str, params))) + '.pkl')
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'cache', folder, f_name)
    return path


def get_data_path(benchmark, params):
    f_name = (benchmark.__class__.__name__[:-6]
              + '_data_' + '_'.join(list(map(str, params))) + '.pkl')
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'cache', 'tmp', f_name)
    return path


def clear_tmp():
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'cache', 'tmp')

    list(map(os.remove, (os.path.join(path, f) for f in os.listdir(path))))


class Benchmark(ABC):
    timer = timeit.default_timer  # wall time
    processes = 1
    timeout = 500

    (profile, n_jobs_vals, save_estimators, save_folder, base_folder,
     bench_predict, bench_transform) = get_from_config()

    if profile == 'fast':
        warmup_time = 0
        repeat = 1
        number = 1
        min_run_count = 1
        data_size = 'small'
    elif profile == 'regular':
        warmup_time = 1
        repeat = (3, 100, 30)
        data_size = 'small'
    elif profile == 'large_scale':
        warmup_time = 1
        repeat = 3
        number = 1
        data_size = 'large'

    @property
    @classmethod
    @abstractmethod
    def params(self):
        pass


class Estimator(ABC):
    @abstractmethod
    def setup_cache_(self, params):
        pass

    def setup_cache(self):
        clear_tmp()

        param_grid = list(itertools.product(*self.params))

        for params in param_grid:
            data, estimator = self.setup_cache_(params) or (None, None)
            if data is None:
                continue

            data_path = get_data_path(self, params)
            with open(data_path, 'wb') as f:
                pickle.dump(data, f)

            X, _, y, _ = data
            estimator.fit(X, y)

            est_path = get_estimator_path(self, Benchmark.save_folder,
                                          params, Benchmark.save_estimators)
            with open(est_path, 'wb') as f:
                pickle.dump(estimator, f)

    def setup(self, *params):
        if hasattr(self, 'setup_'):
            self.setup_(params)

        data_path = get_data_path(self, params)
        with open(data_path, 'rb') as f:
            self.X, self.X_val, self.y, self.y_val = pickle.load(f)

        est_path = get_estimator_path(self, Benchmark.save_folder,
                                      params, Benchmark.save_estimators)
        with open(est_path, 'rb') as f:
            self.estimator = pickle.load(f)

        self.make_scorers()

    def time_fit(self, *args):
        self.estimator.fit(self.X, self.y)

    def peakmem_fit(self, *args):
        self.estimator.fit(self.X, self.y)

    def track_train_score(self, *args):
        if hasattr(self.estimator, 'predict'):
            y_pred = self.estimator.predict(self.X)
        else:
            y_pred = None
        return float(self.train_scorer(self.y, y_pred))

    def track_test_score(self, *args):
        if hasattr(self.estimator, 'predict'):
            y_val_pred = self.estimator.predict(self.X_val)
        else:
            y_val_pred = None
        return float(self.test_scorer(self.y_val, y_val_pred))


class Predictor(ABC):
    if Benchmark.bench_predict:
        def time_predict(self, *args):
            self.estimator.predict(self.X)

        def peakmem_predict(self, *args):
            self.estimator.predict(self.X)

        if Benchmark.base_folder is not None:
            def track_same_prediction(self, *args):
                file_path = get_estimator_path(self, Benchmark.base_folder,
                                               args, True)
                with open(file_path, 'rb') as f:
                    estimator_base = pickle.load(f)

                y_val_pred_base = estimator_base.predict(self.X_val)
                y_val_pred = self.estimator.predict(self.X_val)

                return np.allclose(y_val_pred_base, y_val_pred)

    @property
    @classmethod
    @abstractmethod
    def params(self):
        pass


class Transformer(ABC):
    if Benchmark.bench_transform:
        def time_transform(self, *args):
            self.estimator.transform(self.X)

        def peakmem_transform(self, *args):
            self.estimator.transform(self.X)

        if Benchmark.base_folder is not None:
            def track_same_transform(self, *args):
                file_path = get_estimator_path(self, Benchmark.base_folder,
                                               args, True)
                with open(file_path, 'rb') as f:
                    estimator_base = pickle.load(f)

                X_val_t_base = estimator_base.transform(self.X_val)
                X_val_t = self.estimator.transform(self.X_val)

                return np.allclose(X_val_t_base, X_val_t)

    @property
    @classmethod
    @abstractmethod
    def params(self):
        pass
