"""
===========================================
Partially Observable Markov Process Diagram
===========================================

This script generates the Graphical Model used to represent a Partially
Observable Markov Process.  In this discrete-time model, the unknown hidden
state is represented by :math:`x_t`, and the known observations are represented
by :math:`z_t`.
"""
from matplotlib.patches import Circle, Arrow
from matplotlib import rc
import pylab as pl

boxcolor = 'w'

# make the figure/axes to be drawn on
fig = pl.figure(figsize=(5, 3))
ax = pl.axes((0, 0, 1, 1), xticks=[], yticks=[], frameon=False)
ax.set_xlim(0, 5)
ax.set_ylim(0, 3)

# make the shapes
patches = [
    Circle((1, 2), 0.33, fc=boxcolor),
    Circle((2, 2), 0.33, fc=boxcolor),
    Circle((2.8, 2), 0.05, fc=boxcolor),
    Circle((3.0, 2), 0.05, fc=boxcolor),
    Circle((3.2, 2), 0.05, fc=boxcolor),
    Circle((4, 2), 0.33, fc=boxcolor),
    Arrow(1.4, 2, 0.2, 0, width=0.5, fc=boxcolor),
    Arrow(2.4, 2, 0.2, 0, width=0.5, fc=boxcolor),
    Arrow(3.4, 2, 0.2, 0, width=0.5, fc=boxcolor),
    Circle((1, 1), 0.33, fc=boxcolor),
    Circle((2, 1), 0.33, fc=boxcolor),
    Circle((4, 1), 0.33, fc=boxcolor),
    Arrow(1, 1.6, 0, -0.2, width=0.5, fc=boxcolor),
    Arrow(2, 1.6, 0, -0.2, width=0.5, fc=boxcolor),
    Arrow(4, 1.6, 0, -0.2, width=0.5, fc=boxcolor),
]
for p in patches:
    ax.add_patch(p)

# make the text
rc('text', usetex=True)
pl.text(1, 2, '$x_0$', ha='center', va='center', fontsize=20)
pl.text(2, 2, '$x_1$', ha='center', va='center', fontsize=20)
pl.text(4, 2, '$x_{T-1}$', ha='center', va='center', fontsize=20)
pl.text(1, 1, '$z_0$', ha='center', va='center', fontsize=20)
pl.text(2, 1, '$z_1$', ha='center', va='center', fontsize=20)
pl.text(4, 1, '$z_{T-1}$', ha='center', va='center', fontsize=20)

pl.show()
