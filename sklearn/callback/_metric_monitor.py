# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import inspect
from multiprocessing import Manager

from sklearn.callback._callback_context import get_context_path


class MetricMonitor:
    """Callback that monitors a metric for each iterative steps of an estimator.

    The specified metric function is called on the target values `y` and the predicted
    values on the samples `y_pred = estimator.predict(X)` at each iterative step of the
    estimator.

    Parameters
    ----------
    metric : function
        The metric to compute.
    metric_params : dict or None, default=None
        Additional keyword arguments for the metric function.
    on_validation : bool, default=True
        Whether to compute the metric on validation data (if True) or training data
        (if False).
    """

    def __init__(self, metric, metric_params=None, on_validation=True):
        self.on_validation = on_validation
        self.metric_params = metric_params or dict()
        if metric_params is not None:
            valid_params = inspect.signature(metric).parameters
            invalid_params = [arg for arg in metric_params if arg not in valid_params]
            if invalid_params:
                raise ValueError(
                    f"The parameters '{invalid_params}' cannot be used with the"
                    f" function {metric.__module__}.{metric.__name__}."
                )
        self.metric_func = metric
        self._shared_mem_log = Manager().list()

    def on_fit_begin(self, estimator):
        if not hasattr(estimator, "predict"):
            raise ValueError(
                f"Estimator {estimator.__class__} does not have a predict method, which"
                " is necessary to use a MetricMonitor callback."
            )

    def on_fit_task_end(
        self, estimator, context, data, from_reconstruction_attributes, **kwargs
    ):
        X, y = (
            (data["X_val"], data["y_val"])
            if self.on_validation
            else (data["X_train"], data["y_train"])
        )
        if X is not None and y is not None:
            y_pred = from_reconstruction_attributes().predict(X)
            metric_value = self.metric_func(y, y_pred, **self.metric_params)
        else:
            metric_value = None
        log_item = {self.metric_func.__name__: metric_value}
        for depth, ctx in enumerate(get_context_path(context)):
            if depth == 0:
                timestamp = ctx.init_time.strftime("%Y-%m-%d_%H:%M:%S.%f")
                log_item["_run"] = (
                    f"{ctx.estimator_name}_{ctx.estimator_id}_{timestamp}"
                )
            prev_task_str = (
                f"{ctx.source_estimator_name}_{ctx.source_task_name}|"
                if ctx.source_estimator_name is not None
                else ""
            )
            log_item[f"{depth}_{prev_task_str}{ctx.estimator_name}_{ctx.task_name}"] = (
                ctx.task_id
            )

        self._shared_mem_log.append(log_item)

    def on_fit_end(self, estimator, context):
        pass

    def get_logs(self):
        """Get the logged values.

        A dict of structured dataframes is returned if pandas is installed, otherwise a
        list of dicts is returned.

        If pandas is installed, a dictionary with run ids as keys and mulit-index
        Dataframes with indices corresponding to the task tree as values. The run ids
        are strings of the form : "<estimator name>_<estimator object id>_<timestamp>".

        If pandas is not available, a list of dictionaries. Each dictionary contains :
            "<name of the metric>": <metric value>
            "_run": <estimator name>_<estimator object id>_<timestamp>
        and several pairs of key values describing the task tree :
            "<task depth>_<estimator name>_<task name>": <task id>

        Returns
        -------
        dict of pandas.DataFrame or list of dict
            The logged values.
        """
        logs = list(self._shared_mem_log)

        # If pandas is installed, a structured dataframe is returned, otherwise the raw
        # list of dicts.
        try:
            import pandas as pd

            logs = pd.DataFrame(logs)
            log_dict = {}
            if not logs.empty:
                for run_id in logs["_run"].unique():
                    run_log = logs.loc[logs["_run"] == run_id].copy()
                    # Drop columns that correspond to other runs task_id and are filled
                    # with NaNs, drop the run column, but always keep the metric column.
                    columns_to_keep = ~(run_log.isnull().all())
                    columns_to_keep["_run"] = False
                    columns_to_keep[self.metric_func.__name__] = True
                    run_log = run_log.loc[:, columns_to_keep]
                    log_dict[run_id] = run_log.set_index(
                        [
                            col
                            for col in run_log.columns
                            if col != self.metric_func.__name__
                        ]
                    ).sort_index()
            return log_dict

        except ImportError:
            return logs
