import time
from sklearn.decomposition import FastICA
import numpy as np
import pandas as pd
import streamlit as st
import altair as alt
import argparse

parser = argparse.ArgumentParser(
    description="Determine whether to save/load the dataframe."
)
parser.add_argument(
    "--save",
    type=str,
    default="",
    metavar="s",
    help="Saves the dataframe to the path if provided",
)
parser.add_argument(
    "--load",
    type=str,
    default="",
    metavar="l",
    help="Loads the dataframe at the path if provided",
)
args = parser.parse_args()


def main() -> None:

    df = None
    if args.load:
        df = on_load(args.load)
    else:
        X_shapes = []
        for i in range(10):
            X_shapes.extend(
                (10 * int(10 ** (i / 3)), 10 * int(10 ** (j / 3))) for j in range(10)
            )

        solvers = ("svd", "eigh")
        transformers = {
            s: FastICA(n_components=7, random_state=0, svd_solver=s) for s in solvers
        }
        data = []
        for shape in X_shapes:
            X = np.random.rand(*shape)
            for s in transformers:
                start = time.time()
                transformers[s].fit_transform(X)
                data.append(
                    {"shape": str(shape), "solver": s, "time": time.time() - start}
                )
        df = pd.DataFrame(data)

    if args.save and not args.load:
        df.to_csv(args.save)
        print(f"Dataframe saved to {args.save}")

    chart = (
        alt.Chart(df, width=300)
        .mark_bar()
        .encode(x="solver", y="time", column="shape")
        .properties(title="time by shape")
    ).resolve_scale(y="independent")

    st.altair_chart(chart)

    pt = pd.pivot_table(df, values="time", index=["shape"], columns=["solver"])
    ratio_df = pd.DataFrame(pt["eigh"] / pt["svd"], columns=["eigh/svd time"])
    st.table(ratio_df.sort_values(by="eigh/svd time"))


def on_load(pth: str) -> pd.DataFrame:
    print(f"Dataframe loaded from {args.load}")
    return pd.read_csv(pth)


if __name__ == "__main__":
    main()
