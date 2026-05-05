# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import datetime
import uuid
from collections import defaultdict
from dataclasses import dataclass

from sklearn.callback._callback_context import get_context_path
from sklearn.callback._callback_support import get_callback_manager
from sklearn.utils._optional_dependencies import check_pandas_support
from sklearn.utils._param_validation import StrOptions, validate_params


@dataclass
class ScoringMonitorLog:
    """Log for one run of a scoring monitor.

    The recorded scores are accessed through the `data` attribute, as a list of dicts,
    or the `data_as_pandas` attribute, as a Pandas DataFrame. In the former case, each
    dict corresponds to one row of the corresponding DataFrame and contains
    column_name -> value pairs. The columns are structured as follows:

    - `task_id_path`: tuple containing the task ids from the root task to the
      task for which the score was computed. Each value in this column is
      unique.
    - `parent_task_id_path`: tuple containing the task ids from the root to the
      parent task. It can be used to group scores from tasks that have the same
      parent task.
    - `estimator_name`: the name of the estimator.
    - `task_name`: the name of the task.
    - `task_id`: the id of the task.
    - `sequential_subtasks`: whether the task has sequential subtasks.
    - A column for each score name that was passed as `scoring` parameter.

    Attributes
    ----------
    run_id : uuid.UUID
        The unique identifier for the run.

    estimator_name : str
        The name of the estimator for the run.

    timestamp : datetime.datetime
        The timestamp of the start of the run.

    data : list[dict]
        The recorded scores for the run.

    data_as_pandas : pandas.DataFrame
        The recorded scores for the run as a Pandas DataFrame.
    """

    run_id: uuid.UUID
    estimator_name: str
    timestamp: datetime.datetime
    data: list[dict]

    _data_as_pandas = None

    @property
    def data_as_pandas(self):
        pd = check_pandas_support(f"`{self.__class__.__name__}.data_as_pandas`")
        if self._data_as_pandas is None:
            self._data_as_pandas = pd.DataFrame(self.data)
        return self._data_as_pandas

    def __repr__(self):
        return (
            f"ScoringMonitorLog(run_id={self.run_id}, "
            f"estimator_name={self.estimator_name}, "
            f"timestamp={self.timestamp})"
        )


class ScoringMonitor:
    """Callback that monitors a score for each iterative step of an estimator.

    The specified scorer is called on the training data at each iterative step of the
    estimator, and the score is logged by the callback. The logs can be retrieved
    through the `get_logs` method.

    Parameters
    ----------
    scoring : str, callable, list, tuple, dict or None
        The scoring method to use to monitor the model.

        If `scoring` represents a single score, one can use:

        - a single string (see :ref:`scoring_string_names`);
        - a callable (see :ref:`scoring_callable`) that returns a single value;
        - `None`, the `estimator`'s
          :ref:`default evaluation criterion <scoring_api_overview>` is used.

        If `scoring` represents multiple scores, one can use:

        - a list or tuple of unique strings;
        - a callable returning a dictionary where the keys are the metric
          names and the values are the metric scores;
        - a dictionary with metric names as keys and callables as values.
    """

    @validate_params(
        {"scoring": [str, callable, list, tuple, dict, None]},
        prefer_skip_nested_validation=True,
    )
    def __init__(self, *, scoring):
        self.scoring = scoring
        # Turn the scorer into a MultimetricScorer for convenience
        from sklearn.metrics._scorer import _BaseScorer

        if isinstance(self.scoring, str):
            self.scoring = [self.scoring]
        if callable(self.scoring) and isinstance(self.scoring, _BaseScorer):
            self.scoring = {"score": self.scoring}

        self._shared_log = get_callback_manager().list()
        self._estimator_scorers = {}

    def setup(self, estimator, context):
        # A scorer per estimator is needed to avoid race conditions when the callback is
        # set on different estimators and the scorer is the estimator's default scorer.
        if context.estimator_name not in self._estimator_scorers:
            from sklearn.metrics import check_scoring

            self._estimator_scorers[context.estimator_name] = check_scoring(
                estimator, self.scoring
            )

    def teardown(self, estimator, context):
        pass

    def on_fit_task_begin(self, estimator, context):
        pass

    def on_fit_task_end(
        self,
        estimator,
        context,
        *,
        X=None,
        y=None,
        fitted_estimator=None,
        metadata=None,
    ):
        if fitted_estimator is None:
            return

        context_path = get_context_path(context)

        root_context = context_path[0]
        run_id = root_context.root_uuid
        run_info = {
            "timestamp": root_context.init_time.strftime("UTC%Y-%m-%d-%H:%M:%S.%f"),
            "estimator_name": root_context.estimator_name,
        }

        task_info_path = [
            {
                "estimator_name": ctx.estimator_name,
                "task_name": ctx.task_name,
                "task_id": ctx.task_id,
                "sequential_subtasks": ctx.sequential_subtasks,
            }
            for ctx in context_path
        ]

        scores = {}
        if X is not None and y is not None:
            scorer = self._estimator_scorers[context.estimator_name]
            scores.update(scorer(fitted_estimator, X, y, **metadata))

        self._shared_log.append((run_id, run_info, task_info_path, scores))

    @validate_params(
        {
            "select": [StrOptions({"all", "most_recent"})],
            "as_frame": ["boolean"],
            "include_lineage": ["boolean"],
        },
        prefer_skip_nested_validation=True,
    )
    def get_logs(self, select="most_recent", include_lineage=False):
        """Retrieve the logged scores.

        Log entries are grouped by runs, which are the outermost enclosing fit calls. If
        the estimator this callback is registered on is wrapped in meta-estimators, a
        run corresponds to one fit of the outermost meta-estimator. If it is not wrapped
        in a meta-estimator, a run simply corresponds to a single fit of the estimator.

        For a given run, the scores are logged in a :class:`ScoringMonitorLog` object
        containing:

        - `run_id`: a unique identifier for the run;
        - `estimator_name`: the name of the (meta-)estimator of the run;
        - `timestamp`: the timestamp of the start of the run;
        - `data`: the recorded scores for the run. Each score value is associated
          with the context of the task for which the score was computed;
        - `data_as_pandas`: the recorded scores as a Pandas DataFrame.

        See :class:`ScoringMonitorLog` for more details about the structure of the
        recorded scores.

        Parameters
        ----------
        select : {"all", "most_recent"}, default="most_recent"
            Which log run to return:

            - `"all"`: return the logged scores for all runs;
            - `"most_recent"`: return the logged scores for the most recent run.

        include_lineage : bool, default=False
            Whether to include lineage information of the tasks in the log.

            If set to True, the log contains extra rows for each task that is an
            ancestor of a task for which the score was computed. These extra rows can
            be used to retrieve the context of all ancestor tasks of a given task for
            which the score was computed. For these extra rows, there are no score
            entries if `as_frame` is False, or NaN values if `as_frame` is True.

        Returns
        -------
        logs : :class:`ScoringMonitorLog` or list of :class:`ScoringMonitorLog`
            The logged scores. If `select=="most_recent"`, returns a single
            :class:`ScoringMonitorLog` object. Otherwise, returns the list of all
            run logs.
        """
        logs = defaultdict(lambda: {"data": []})
        run_to_task_id_path = defaultdict(set)

        if len(shared_log := self._shared_log) == 0:
            raise ValueError(
                "No logs to retrieve. No scores were computed during the runs or the "
                "estimator is not fitted yet"
            )

        # group logs by run
        for run_id, run_info, task_info_path, scores in shared_log:
            logs[run_id].update(run_info)

            task_id_path = tuple(task_info["task_id"] for task_info in task_info_path)

            logs[run_id]["data"].append(
                {
                    "task_id_path": task_id_path,
                    "parent_task_id_path": task_id_path[:-1],
                    "parent_task_info_path": task_info_path[:-1],
                    **task_info_path[-1],
                    **scores,
                }
            )
            run_to_task_id_path[run_id].add(task_id_path)

        for run_id, log in logs.items():
            if include_lineage:
                extra_rows = []
                for row in log["data"]:
                    for i in range(len(row["parent_task_info_path"])):
                        task_info_path = row["parent_task_info_path"][: i + 1]
                        task_id_path = tuple(
                            task_info["task_id"] for task_info in task_info_path
                        )
                        if task_id_path not in run_to_task_id_path[run_id]:
                            extra_rows.append(
                                {
                                    "task_id_path": task_id_path,
                                    "parent_task_id_path": task_id_path[:-1],
                                    "parent_task_info_path": task_info_path[:-1],
                                    **task_info_path[-1],
                                }
                            )
                            run_to_task_id_path[run_id].add(task_id_path)
                log["data"] += extra_rows

            # sort rows by recursive task ids so that tasks of a same parent are grouped
            sorting_key = lambda x: (len(x["task_id_path"]), x["task_id_path"])
            log["data"] = sorted(log["data"], key=sorting_key)

            for row in log["data"]:
                row.pop("parent_task_info_path", None)

        # sort logs by run timestamp and estimator name
        logs = [ScoringMonitorLog(run_id=run_id, **log) for run_id, log in logs.items()]
        logs.sort(key=lambda log: (log.timestamp, log.estimator_name))

        if select == "most_recent":
            return logs[-1]

        return logs
