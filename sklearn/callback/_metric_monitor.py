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
    on : str, default="train_set"
        Which data to compue the metric on. Possible values are "train_set",
        "validation_set" and "both". "train_set" corresponds to using the X and y
        arguments of the fit function, "validation_set" corresponds to using the X_val
        and y_val arguments, "both" corresponds to using both.
    """

    def __init__(self, metric, metric_params=None, on="train_set"):
        # TODO: use a scorer for the metric
        possible_on_values = ("train_set", "validation_set", "both")
        if on not in possible_on_values:
            raise ValueError(
                f"'on' should be one of {possible_on_values}, got {on} instead."
            )
        self.on = on
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
        self._shared_log = Manager().list()

    def on_fit_begin(self, estimator):
        if not hasattr(estimator, "predict"):
            raise ValueError(
                f"Estimator {estimator.__class__} does not have a predict method, which"
                " is necessary to use a MetricMonitor callback."
            )

    def on_fit_task_end(
        self, estimator, context, data, from_reconstruction_attributes, **kwargs
    ):
        # TODO: add a task_info dict in the logs
        reconstructed_est = from_reconstruction_attributes()
        context_path = get_context_path(context)
        if self.on == "train_set" or self.on == "both":
            X, y = None, None
            if "X_train" in data and "y_train" in data:
                X, y = data["X_train"], data["y_train"]
            self._log_item(X, y, "train_set", reconstructed_est, context_path)
        if self.on == "validation_set" or self.on == "both":
            X, y = None, None
            if "X_val" in data and "y_val" in data:
                X, y = data["X_val"], data["y_val"]
            self._log_item(X, y, "validation_set", reconstructed_est, context_path)

    def _log_item(self, X, y, on, reconstructed_est, context_path):
        if X is not None and y is not None:
            y_pred = reconstructed_est.predict(X)
            metric_value = self.metric_func(y, y_pred, **self.metric_params)
        else:
            metric_value = None
        log_item = {self.metric_func.__name__: metric_value, "on": on}
        for depth, ctx in enumerate(context_path):
            if depth == 0:
                timestamp = ctx.init_time.strftime("%Y-%m-%d_%H:%M:%S.%f")
                # TODO: use a UUID in the run id instead or in addition of the timestamp
                log_item["run"] = f"{ctx.estimator_name}_{ctx.estimator_id}_{timestamp}"
            prev_task_str = (
                f"{ctx.source_estimator_name}_{ctx.source_task_name}|"
                if ctx.source_estimator_name is not None
                else ""
            )
            log_item[f"{depth}_{prev_task_str}{ctx.estimator_name}_{ctx.task_name}"] = (
                ctx.task_id
            )

        self._shared_log.append(log_item)

    def on_fit_end(self, estimator, context):
        pass

    def get_logs(self):
        """Get the logged values.

        A list of tuples containing structured dataframes is returned if pandas is
        installed, otherwise a list of dicts of lists is returned. In both cases, if
        there is only one run in the logs, the item is returned directly instead of a
        singleton list.

        If pandas is installed, a list of tuples is returned. Each tuple contains a
        run id and a mulit-index Dataframe with indices corresponding to the task tree
        as values. The run ids are strings of the form : "<estimator name>_<estimator
        object id>_<timestamp>".

        If pandas is not available, a list of dictionaries of lists. Each dictionary
        corresponds to a run and contains :
            "<name of the metric>": list of metric values
            "run": <estimator name>_<estimator object id>_<timestamp>
        and several pairs of key values describing the task tree :
            "<task depth>_<estimator name>_<task name>": list of task id

        Returns
        -------
        dict of pandas.DataFrame or list of dicts of lists
            The logged values.
        """
        logs = list(self._shared_log)

        try:
            import pandas as pd

            logs = pd.DataFrame(logs)
            log_list = []
            if not logs.empty:
                for run_id in logs["run"].unique():
                    run_log = logs.loc[logs["run"] == run_id].copy()
                    # Drop columns that correspond to other runs task_id which are
                    # filled with NaNs.
                    columns_to_keep = ~(run_log.isnull().all())
                    columns_to_keep["run"] = False
                    columns_to_keep[self.metric_func.__name__] = True
                    run_log = run_log.loc[:, columns_to_keep]
                    log_list.append(
                        (
                            run_id,
                            run_log.set_index(
                                [
                                    col
                                    for col in run_log.columns
                                    if col not in (self.metric_func.__name__, "on")
                                ]
                            ).sort_index(),
                        )
                    )
            return log_list[0] if len(log_list) == 1 else log_list

        except ImportError:
            runs_dict_list = []
            for log_item in logs:
                run_name = log_item["run"]
                for run_dict in runs_dict_list:
                    if run_dict["run"] == run_name:
                        for key, val in log_item.items():
                            if key != "run":
                                run_dict[key].append(val)
                        break
                else:
                    run_dict = {"run": run_name}
                    for key, val in log_item.items():
                        if key != "run":
                            run_dict[key] = [val]
                    runs_dict_list.append(run_dict)

            return runs_dict_list[0] if len(runs_dict_list) == 1 else runs_dict_list
