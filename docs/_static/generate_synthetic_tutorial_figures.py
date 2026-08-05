"""Generate fixed-seed, non-human tutorial figures for the beginner guide."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SEED = 20260804
OUT_DIR = Path(__file__).resolve().parent
COLORS = {"Control": "#0072B2", "Case": "#D55E00"}


def _style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#d8dde2", linewidth=0.7, alpha=0.7)


def _smooth_histogram(values, bins, bandwidth=7):
    density, edges = np.histogram(values, bins=bins, density=True)
    x = (edges[:-1] + edges[1:]) / 2
    kernel_x = np.arange(-3 * bandwidth, 3 * bandwidth + 1)
    kernel = np.exp(-0.5 * (kernel_x / bandwidth) ** 2)
    kernel /= kernel.sum()
    return x, np.convolve(density, kernel, mode="same")


def methylation_distribution(rng):
    fig, ax = plt.subplots(figsize=(6.4, 4.1), constrained_layout=True)
    for group, shift in (("Control", 0.00), ("Case", 0.02)):
        low = rng.beta(1.2, 7.0, 28_000)
        high = rng.beta(6.5, 1.5, 22_000)
        values = np.clip(np.concatenate([low, high]) + shift, 0, 1)
        x, density = _smooth_histogram(values, np.linspace(0, 1, 151), bandwidth=4)
        ax.plot(x, density, color=COLORS[group], linewidth=2.2, label=group)
    ax.set(xlabel="CpG methylation beta value", ylabel="Density",
           title="Synthetic methylation distribution")
    ax.set_xlim(0, 1)
    ax.legend(frameon=False)
    _style_axis(ax)
    fig.savefig(OUT_DIR / "tutorial_methylation_distribution.png", dpi=180,
                facecolor="white")
    plt.close(fig)


def fragment_length_distribution(rng):
    fig, ax = plt.subplots(figsize=(6.4, 4.1), constrained_layout=True)
    for group, shift in (("Control", 0.0), ("Case", 2.0)):
        mono = rng.normal(167 + shift, 16, 48_000)
        di = rng.normal(325 + shift, 25, 5_000)
        values = np.clip(np.concatenate([mono, di]), 30, 500)
        x, density = _smooth_histogram(values, np.arange(30, 501, 2), bandwidth=3)
        ax.plot(x, 100 * density, color=COLORS[group], linewidth=2.2, label=group)
    ax.axvline(167, color="#555555", linewidth=1, linestyle="--", label="167 bp")
    ax.set(xlabel="Fragment length (bp)", ylabel="Density (scaled %)",
           title="Synthetic cfDNA fragment-length distribution")
    ax.set_xlim(50, 400)
    ax.legend(frameon=False)
    _style_axis(ax)
    fig.savefig(OUT_DIR / "tutorial_fragment_length_distribution.png", dpi=180,
                facecolor="white")
    plt.close(fig)


def main():
    rng = np.random.default_rng(SEED)
    methylation_distribution(rng)
    fragment_length_distribution(rng)


if __name__ == "__main__":
    main()
