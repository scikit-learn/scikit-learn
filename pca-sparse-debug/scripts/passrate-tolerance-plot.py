import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

csv = sys.argv[1]
plot = sys.argv[2]

df = pd.read_csv(csv)
df['pass'] = df.outcome.apply(lambda x: True if x=='PASSED' else False)
df = df[df.solver != 'auto']

def passrate_by(x):
    passes = df.groupby(x)['pass']
    counts = passes.count()
    sums = passes.sum()
    return sums / counts

gb = passrate_by(['solver', 'rtol'])
def formatter(x):
    if x == 1:
        return '1e-00'
    else:
        return format(x, ".0e")
gb.index = gb.index.set_levels(map(formatter, gb.index.levels[1]), level=1)
fig, axes = plt.subplots(1, 4, figsize=(24, 4), dpi=300)
fig.suptitle('Test pass rate by solver and relative tolerance', fontsize=20)
for solver, ax in zip(gb.index.levels[0], axes.flat):
    bar = gb[solver].plot.bar(ax=ax)
    bar.set_title(solver, pad=15, fontsize=20)
    bar.set_xlabel('relative tolerance', fontsize=14)
    bar.set_ylabel('pass rate', fontsize=14)
    bar.set_ylim(bottom=0, top=1)
    bar.tick_params(axis='both', which='major', labelsize=14)
    ax.bar_label(ax.containers[0], fmt="%.2g")

fig.tight_layout()
fig.savefig(plot, facecolor='white', transparent=False)
