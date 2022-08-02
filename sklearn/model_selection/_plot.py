import numpy as np

from . import learning_curve
from ..utils import check_matplotlib_support


class LearningCurveDisplay:
    def __init__(self, *, train_sizes, train_scores, test_scores, score_name):
        self.train_sizes = train_sizes
        self.train_scores = train_scores
        self.test_scores = test_scores
        self.score_name = score_name

    def plot(
        self,
        ax=None,
        *,
        is_score=True,
        score_name=None,
        log_scale=False,
        std_display_style="errorbar",
        line_kw=None,
        fill_between_kw=None,
        errorbar_kw=None,
    ):
        check_matplotlib_support(f"{self.__class__.__name__}.plot")

        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()

        if is_score:
            train_scores, test_scores = self.train_scores, self.test_scores
            label = "Score"
        else:
            train_scores, test_scores = -self.train_scores, -self.test_scores
            label = "Error"

        if std_display_style == "errorbar":
            if errorbar_kw is None:
                errorbar_kw = {}
            errorbar_train = ax.errorbar(
                x=self.train_sizes,
                y=train_scores.mean(axis=1),
                yerr=train_scores.std(axis=1),
                label=f"Training {label}",
                **errorbar_kw,
            )
            errorbar_test = ax.errorbar(
                x=self.train_sizes,
                y=test_scores.mean(axis=1),
                yerr=test_scores.std(axis=1),
                label=f"Testing {label}",
                **errorbar_kw,
            )
            self.errorbar_ = [errorbar_train, errorbar_test]
            self.line_, self.fill_between_ = None, None
        elif std_display_style == "fill_between":
            if line_kw is None:
                line_kw = {}
            if fill_between_kw is None:
                fill_between_kw = {}
            default_fill_between_kw = {"alpha": 0.5}
            fill_between_kw = {**default_fill_between_kw, **fill_between_kw}
            line_train = ax.plot(
                self.train_sizes,
                train_scores.mean(axis=1),
                label=f"Training {label}",
                **line_kw,
            )
            fill_between_train = ax.fill_between(
                x=self.train_sizes,
                y1=train_scores.mean(axis=1) - train_scores.std(axis=1),
                y2=train_scores.mean(axis=1) + train_scores.std(axis=1),
                **fill_between_kw,
            )
            line_test = ax.plot(
                self.train_sizes,
                test_scores.mean(axis=1),
                label=f"Testing {label}",
                **line_kw,
            )
            fill_between_test = ax.fill_between(
                x=self.train_sizes,
                y1=test_scores.mean(axis=1) - test_scores.std(axis=1),
                y2=test_scores.mean(axis=1) + test_scores.std(axis=1),
                **fill_between_kw,
            )
            self.line_ = line_train + line_test
            self.fill_between_ = [fill_between_train, fill_between_test]
            self.errorbar_ = None
        else:
            raise ValueError(
                f"Unknown std_display_style: {std_display_style}. Should be one of"
                " 'errorbar' or 'fill_between'"
            )

        score_name = self.score_name if score_name is None else score_name

        ax.legend()
        if log_scale:
            ax.set_xscale("log")
        ax.set_xlabel("Number of samples in the training set")
        ax.set_ylabel(f"{score_name}")

        self.ax_ = ax
        self.figure_ = ax.figure
        return self

    @classmethod
    def from_estimator(
        cls,
        estimator,
        X,
        y,
        *,
        groups=None,
        train_sizes=np.linspace(0.1, 1.0, 5),
        cv=None,
        scoring=None,
        exploit_incremental_learning=False,
        n_jobs=None,
        pre_dispatch="all",
        verbose=0,
        shuffle=False,
        random_state=None,
        error_score=np.nan,
        fit_params=None,
        ax=None,
        is_score=True,
        score_name=None,
        log_scale=False,
        std_display_style="errorbar",
        line_kw=None,
        fill_between_kw=None,
        errorbar_kw=None,
    ):
        check_matplotlib_support(f"{cls.__name__}.from_estimator")

        score_name = "Score" if score_name is None else score_name

        train_sizes, train_scores, test_scores = learning_curve(
            estimator,
            X,
            y,
            groups=groups,
            train_sizes=train_sizes,
            cv=cv,
            scoring=scoring,
            exploit_incremental_learning=exploit_incremental_learning,
            n_jobs=n_jobs,
            pre_dispatch=pre_dispatch,
            verbose=verbose,
            shuffle=shuffle,
            random_state=random_state,
            error_score=error_score,
            return_times=False,
            fit_params=fit_params,
        )

        viz = cls(
            train_sizes=train_sizes,
            train_scores=train_scores,
            test_scores=test_scores,
            score_name=score_name,
        )
        return viz.plot(
            ax=ax,
            is_score=is_score,
            log_scale=log_scale,
            std_display_style=std_display_style,
            line_kw=line_kw,
            fill_between_kw=fill_between_kw,
            errorbar_kw=errorbar_kw,
        )
