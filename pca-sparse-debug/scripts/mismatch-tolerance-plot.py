import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

csv = sys.argv[1]
plot = sys.argv[2]

df = pd.read_csv(csv)
df = df[df.solver != 'auto']

df['rate'] = df.bad/df.total

def mismatch_by(x):
    gb = df.groupby(x)
    return gb['rate'].mean()

gb = mismatch_by(['solver', 'rtol'])
gb = gb.reindex(pd.MultiIndex.from_product(gb.index.levels)).fillna(0)

def formatter(x):
    if x == 1:
        return '1e-00'
    else:
        return format(x, ".0e")
gb.index = gb.index.set_levels(map(formatter, gb.index.levels[1]), level=1)
solvers = gb.index.levels[0]
fig, axes = plt.subplots(1, len(solvers), figsize=(24, 4), dpi=300)
fig.suptitle('Elementwise mismatch rate by solver and relative tolerance', fontsize=20)
for solver, ax in zip(solvers, axes.flat):
    bar = gb[solver].plot.bar(ax=ax)
    bar.set_title(solver, pad=15, fontsize=20)
    bar.set_xlabel('relative tolerance', fontsize=14)
    bar.set_ylabel('mismatch rate', fontsize=14)
    bar.set_ylim(bottom=0, top=1)
    bar.tick_params(axis='both', which='major', labelsize=14)
    ax.bar_label(ax.containers[0], fmt="%.2g")

fig.tight_layout()
fig.savefig(plot, facecolor='white', transparent=False)
