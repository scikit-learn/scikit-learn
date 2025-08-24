import numpy as np
import time
import pandas as pd
from sklearn.linear_model import QuantileRegressor
from sklearn.metrics import mean_pinball_loss
from tqdm import tqdm
import warnings
import math
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid")


class QuantileRegressorBenchmark(QuantileRegressor):
    def __init__(self, N, k, rng, solver, solver_options=None, n_iter=5):
        super().__init__(solver=solver, solver_options=solver_options, alpha=0)

        self.N = N
        self.k = k
        self.n_iter = n_iter
        self.solver = solver

        self.tol = min(solver_options.values())

        if solver_options is not None:
            solver_options["maxiter"] = solver_options.get("max_iter", None)

        self.rng = rng
        self.y, self.X = self.get_data(N=self.N, k=self.k)

        self.start_time = time.time()
        for i in range(self.n_iter):
            self.fit(X=self.X, y=self.y)

        self.end_time = time.time()
        self.get_pinball_loss()
        self.collect()

    def get_data(self, N, k):
        X = self.rng.normal(size=(N, k))
        y = 1 + 2 * X[:, 0] + 3 * X[:, 1] + self.rng.normal(0, 1, N)
        return y, X

    def get_pinball_loss(self):
        y = self.y
        yhat = self.predict(X=self.X)
        return mean_pinball_loss(y, yhat, alpha=self.quantile)

    def collect(self):
        return {
            "solver": self.solver,
            "N": self.N,
            "k": self.k,
            "coef": self.coef_[0],
            "nit": self.n_iter_,
            "pinball_loss": self.get_pinball_loss(),
            "time": (self.end_time - self.start_time) / self.n_iter,
            "estimation_n_iter": self.n_iter,
            "tol": self.tol,
        }


def _run_qr(solver, N, k, iter_cap, seed=1):
    rng = np.random.default_rng(seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return QuantileRegressorBenchmark(
            solver=solver,
            N=N,
            k=k,
            rng=rng,
            n_iter=1,
            # SciPy's linprog expects 'maxiter' (ensure it's an int)
            solver_options={"maxiter": int(iter_cap)},
        ).collect()


def _run_benchmark(N_list, k_list):
    # first staircase part

    head = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

    # geometric tail: keep halving until 1e-10
    tail = []
    val = 0.5
    while val >= 1e-10:
        tail.append(val)
        val *= 0.5
    tail = np.array(tail)

    head = [10, 9, 8, 7, 6]

    tols = np.concatenate([head, tail])

    all_list = []

    for N in N_list:
        for k in k_list:
            simplex_converged = False
            fn_converged = False

            for tol in tqdm(tols):
                with warnings.catch_warnings():
                    warnings.simplefilter(
                        "ignore"
                    )  # ignore *all* warnings inside the block
                    qr_highs = QuantileRegressorBenchmark(
                        solver="highs",
                        N=N,
                        k=k,
                        rng=np.random.default_rng(1),
                        n_iter=5,
                        solver_options={
                            "ipm_optimality_tolerance": tol,
                            "primal_feasibility_tolerance": tol,
                            "dual_feasibility_tolerance": tol,
                        },
                    ).collect()
                    all_list.append(qr_highs)

            for iterations in range(1, 50):
                with warnings.catch_warnings():
                    warnings.simplefilter(
                        "ignore"
                    )  # ignore *all* warnings inside the block

                    qr_fn = QuantileRegressorBenchmark(
                        solver="frisch-newton",
                        N=N,
                        k=k,
                        rng=np.random.default_rng(1),
                        n_iter=5,
                        solver_options={"max_iter": iterations},
                    ).collect()

                all_list.append(qr_fn)

    res = pd.DataFrame(all_list)
    res.to_csv("fn_vs_highs_benchmark.csv")

    return res


def plot_convergence(df, xcol, xlabel, title, fname, logx=True):
    df = df.copy()
    df["is_highs"] = df["solver"].str.contains("highs")

    g = sns.FacetGrid(
        df,
        col="N",
        row="k",
        sharex=True,
        sharey=True,
        margin_titles=True,
        height=3.2,
        aspect=1.2,
    )

    # Plot line traces
    g.map_dataframe(
        sns.lineplot,
        x=xcol,
        y="pinball_loss",
        hue="solver",
        style="solver",
        markers=True,
        dashes=False,
        estimator=None,
        legend="full",   # let seaborn create legend handles
    )

    # Log-scale for x-axis if requested
    if logx:
        for ax in g.axes.flat:
            ax.set_xscale("log")

    g.set_axis_labels(xlabel, "Mean pinball loss")
    g.set_titles(col_template="N={col_name}", row_template="k={row_name}")
    g.fig.suptitle(title, y=1.03, fontsize=14)

    # Place single legend to the right of the plots
    g.add_legend(title="Solver")  # uses hue="solver"
    g._legend.set_bbox_to_anchor((1.02, 0.5))  # shift outside right
    g._legend.set_frame_on(False)

    g.savefig(fname, dpi=300, bbox_inches="tight")
    plt.close(g.fig)
    print(f"Saved {fname}")


if True:
    N_list = [1000, 5000]
    k_list = [10, 20]
    df = _run_benchmark(N_list=N_list, k_list=k_list)
    # --- Usage ---
    plot_convergence(
        df,
        xcol="nit",
        xlabel="Iterations (log scale)",
        title="Convergence: Pinball loss vs iterations",
        fname=f"pinball_vs_iterations_{N_list[0]}.png",
        logx=True,
    )

    plot_convergence(
        df,
        xcol="time",
        xlabel="Wall-clock time [s] (log scale)",
        title="Convergence: Pinball loss vs time",
        fname=f"pinball_vs_time_{N_list[0]}.png",
        logx=True,
    )

