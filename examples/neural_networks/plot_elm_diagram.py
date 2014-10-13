"""
================
Plot ELM diagram
================

This example plots a single-hidden layer feedforward network. It also allows
users to generate networks of different number of layers and
neurons by changing the two variables 'symbols' and 'n_neurons'.
"""
print(__doc__)

# Author: Issam H. Laradji <issam.laradji@gmail.com>
# License: BSD 3 clause

import pylab as pl
from matplotlib.patches import Circle, Arrow


def create_layer(ax, symbol, n_neurons, bias, x, y):
    patches = []
    if bias is True:
        patches.append(Circle((x, y + 3 * radius), radius, zorder=1,
                              fc='#CCCCCC'))
        pl.text(x, y + 3 * radius, "$+1$", ha='center', va='center',
                fontsize=fontsize)

    for i in range(n_neurons):
        patches.append(Circle((x, y - i * (3 * radius)), radius,
                              zorder=1, fc='#CCCCCC'))

        neuron_symbol = symbol
        if n_neurons != 1:
            neuron_symbol += "$_" + str(i + 1) + "$"

        pl.text(x, y - i * (3 * radius), neuron_symbol,
                ha='center', va='center', fontsize=fontsize)

    for p in patches:
        ax.add_patch(p)

    return patches


def create_arrows(ax, prev_patches, curr_patches):
    for prev in prev_patches:
        for curr in curr_patches:
            dx = curr.center[0] - prev.center[0] - 2 * radius
            dy = curr.center[1] - prev.center[1]
            ax.add_patch(Arrow(prev.center[0] + prev.radius, prev.center[1],
                               dx, dy, antialiased=True, fc='#88CCFF',
                               width=0.05))

# Change 'symbols' and 'n_neurons' for generating different networks
symbols = ["$x$", "$h$", "$f(x)$"]
n_neurons = [3, 3, 1]

assert len(symbols) == len(n_neurons)

radius = 0.6
fontsize = 50 * radius

n_layers = len(symbols)
x_size = n_layers * radius * 3.5
y_size = n_neurons[0] * radius * 4

fig = pl.figure(figsize=(x_size, y_size), facecolor='w')
ax = pl.axes((0, 0, 1, 1), xticks=[], yticks=[], frameon=False)
ax.set_xlim(0, x_size)
ax.set_ylim(0, y_size)

start_x = radius * 2
start_y = y_size - radius * 4.5

rows = start_y + 3 * radius
prev_patches = create_layer(ax, symbols[0], n_neurons[0], True, start_x,
                            start_y)
for i in range(1, n_layers):
    start_y = rows - (rows - (n_neurons[i] * 2 * radius)) / 2
    curr_patches = create_layer(ax, symbols[i], n_neurons[i], False,
                                start_x + 3 * radius, start_y)
    create_arrows(ax, prev_patches, curr_patches)

    prev_patches = curr_patches
    start_x += 3 * radius

pl.show()
