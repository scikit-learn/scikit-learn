# Code to produce the plot poisson_gamma_tweedie_distributions.png

import matplotlib.pyplot as plt
import numpy as np
import tweedie  # see https://pypi.org/project/tweedie/
from scipy import stats

fig, ax = plt.subplots(1, 3, figsize=(12, 3.5), constrained_layout=True)
x = np.arange(0, 11)
mu_grid = [0.5, 1, 3, 5]

# Poisson
for mu in mu_grid:
    ax[0].plot(x, stats.poisson.pmf(x, mu=mu), marker=".", label=rf"$\mu$={mu:.1f}")
ax[0].set_title("Poisson distribution")
ax[0].set_ylim(0, 1)

x = np.linspace(0.01, 10, 1000)

# Tweedie
for mu in mu_grid:
    color = next(ax[1]._get_lines.prop_cycler)["color"]
    ax[1].plot(
        x,
        tweedie.tweedie.pdf(x, mu=mu, p=1.5, phi=1.0),
        label=rf"$\mu$={mu:.1f}",
        color=color,
    )
    # point measure at 0 for 1<=p<2
    ax[1].plot(
        0,
        tweedie.tweedie.pdf(0, mu=mu, p=1.5, phi=1.0),
        marker=".",
        markersize=8,
        color=color,
    )
    ax[1].set_ylim(0, 1)
ax[1].set_title("Tweedie distribution (p=1.5)")

# Gamma
for k in mu_grid:
    ax[2].plot(x, stats.gamma.pdf(x, k, scale=1), label=rf"$\mu$={k:.1f}")
    ax[2].set_ylim(0, 1)
ax[2].set_title("Gamma distribution")

for axi in ax:
    axi.set_xlabel("Y")
    axi.set_xlim(-0.5, 10)
    axi.grid(which="both")
    axi.set_ylabel("PDF(Y)")
    axi.set_facecolor("white")

ax[1].legend()
fig.savefig("poisson_gamma_tweedie_distributions.png", facecolor="None")
