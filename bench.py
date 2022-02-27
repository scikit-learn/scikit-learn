import time
import numpy as np
import pandas as pd
import argparse
from scipy import linalg

# import streamlit as st
# import altair as alt

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
        for i in range(3):
            X_shapes.extend((int(10 ** (2 + i)), int(10 ** (1 + j))) for j in range(3))

        solvers = {
            "svd": linalg.svd,
            "eigh": linalg.eigh,
        }
        total_reps = len(solvers) * len(X_shapes)
        count = 0
        data = []
        for shape in X_shapes:
            XT = np.random.rand(*shape).T

            count += 1
            start = time.time()
            solvers["svd"](XT)
            print(f"Progress: {count}/{total_reps}")
            svd_time = time.time() - start

            count += 1
            start = time.time()
            solvers["eigh"](XT.dot(XT.T))
            print(f"Progress: {count}/{total_reps}")
            eigh_time = time.time() - start

            data.append({"shape": str(shape), "svd": svd_time, "eigh": eigh_time})
        df = pd.DataFrame(data)

    if args.save and not args.load:
        df.to_csv(args.save, index=False)
        print(f"Dataframe saved to {args.save}")

    """
    chart = (
        alt.Chart(df, width=300)
        .mark_bar()
        .encode(x="shape", y=["svd","eigh"], column="shape")
        .properties(title="time by shape")
    ).resolve_scale(y="independent")

    st.altair_chart(chart)
    """


def on_load(pth: str) -> pd.DataFrame:
    print(f"Dataframe loaded from {args.load}")
    df = pd.read_csv(pth)
    df["shape"] = df["shape"].astype("string")
    return df


if __name__ == "__main__":
    main()
