import sys
import pandas as pd
import matplotlib.pyplot as plt

csv = sys.argv[1]
plot = sys.argv[2]

df = pd.read_csv(csv)
df = df[df.solver != 'auto']

df['rate'] = df.bad/df.total

def mismatch_by(x):
    gb = df.groupby(x)
    return gb['rate'].mean()

fig, axes = plt.subplots(2, 2, figsize=(10, 10), dpi=300)

ax=axes[0][0]
seed = mismatch_by('seed').hist(ax=ax)
seed.set_title('mismatch rate by seed (histogram)')
seed.set_xlabel('mismatch rate')
seed.set_ylabel('seed count')
seed.set_ylim(top=100)
seed.set_xlim(right=1)

ax=axes[0][1]
solver = mismatch_by('solver').plot.bar(ax=ax)
solver.set_title('mismatch rate by solver')
solver.set_xlabel('solver')
solver.set_ylabel('mismatch rate')
ax.bar_label(ax.containers[0], fmt="%.3f")

ax=axes[1][0]
density = mismatch_by('density').plot.bar(ax=ax)
density.set_title('mismatch rate by density')
density.set_xlabel('density')
density.set_ylabel('mismatch rate')
ax.bar_label(ax.containers[0], fmt="%.3f")

ax=axes[1][1]
ncomp = mismatch_by('k').plot.bar(ax=ax)
ncomp.set_title('mismatch rate by number of components')
ncomp.set_xlabel('# components')
ncomp.set_ylabel('mismatch rate')
ax.bar_label(ax.containers[0], fmt="%.3f")

for bp in (solver, density, ncomp):
    bp.set_xticklabels(bp.get_xticklabels(), rotation=0)
    bp.set_ylim(top=1)
fig.tight_layout()
fig.savefig(plot, facecolor='white', transparent=False)
