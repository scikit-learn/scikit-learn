# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections import defaultdict

from sklearn.callback._callback_context import get_context_path
from sklearn.callback._callback_support import get_callback_manager
from sklearn.metrics import check_scoring
from sklearn.metrics._scorer import _BaseScorer
from sklearn.utils._optional_dependencies import check_pandas_support
from sklearn.utils._param_validation import StrOptions, validate_params


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

        context_levels = {"task_path": tuple()}
        for depth, ctx in enumerate(context_path):
            context_levels["task_path"] += (ctx.task_id,)
            for k in ("estimator_name", "task_name", "task_id", "sequential_subtasks"):
                context_levels[f"{k}_depth_{depth}"] = getattr(ctx, k)

        scores = {}
        if X is not None and y is not None:
            scorer = self._estimator_scorers[context.estimator_name]
            scores.update(scorer(fitted_estimator, X, y, **metadata))

        self._shared_log.append((run_id, run_info, context_levels, scores))

    @validate_params(
        {
            "select": [StrOptions({"all", "most_recent"})],
            "as_frame": ["boolean"],
        },
        prefer_skip_nested_validation=True,
    )
    def get_logs(self, select="most_recent", as_frame=False):
        """Get the logged scores.

        The callback records scores for the estimator it is registered to. Log entries
        are grouped by the outermost enclosing fit call (the fit of the outermost
        meta-estimator that is a parent of the estimator the callback is registered to).
        Thus one run corresponds to one top-level fit. If the estimator is not wrapped
        in a meta-estimator, a run simply corresponds to a single fit of the estimator.

        For a given run, the scores are logged in a dict containing:
            - "run_id": a unique identifier for the run;
            - "estimator_name": the name of the estimator of the run;
            - "timestamp": the timestamp of the start of the run;
            - "data": the recorded scores for the run. Each score value is associated
              with the detailed context of the score computation.

        Depending on `as_frame`, "data" is either a Pandas DataFrame or a list of dicts.
        In the latter case, each dict corresponds to one row of the corresponding
        DataFrame and contains column_name -> value pairs. The columns are structured as
        follows:
        - "task_path": tuple containing the task ids of the task and its ancestors for
            which the score was computed;
        - A group of 4 columns for each ancestor task of the task for which the score
            was computed (including itself): "estimator_name_depth_<depth>",
            "task_name_depth_<depth>", "task_id_depth_<depth>",
            "sequential_subtasks_depth_<depth>".
        - A column for each score name the was passed as `scoring` parameter.

        Parameters
        ----------
        select : {"all", "most_recent"}, default="most_recent"
            Which log run to return:

            - `"all"`: return the logged scores for all runs;
            - `"most_recent"`: return the logged scores for the most recent run.

        as_frame : bool, default=False
            Whether to have the individual run logs formatted as Pandas DataFrames. If
            set to False the individual run logs are formatted as lists of dictionaries
            instead.

        Returns
        -------
        logs : dict or list of dict
            The logged scores.
        """
        logs = defaultdict(lambda: {"data": []})

        # group logs by run
        for run_id, run_info, context_levels, scores in list(self._shared_log):
            logs[run_id].update(run_info)
            logs[run_id]["data"].append({**context_levels, **scores})

        # sort logs by run timestamp and estimator name
        logs = [{"run_id": run_id, **log} for run_id, log in logs.items()]
        logs.sort(key=lambda log: (log["timestamp"], log["estimator_name"]))

        if as_frame:
            pd = check_pandas_support(f"`{self.__class__.__name__}.get_logs`")
            for log in logs:
                log["data"] = pd.DataFrame(log["data"])

        if select == "most_recent":
            return logs[-1] if logs else {}

        return logs
