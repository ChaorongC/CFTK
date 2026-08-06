"""Generate fixed-seed, non-human tutorial figures for CFTK documentation."""

from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SEED = 20260804
OUT_DIR = Path(__file__).resolve().parent
DPI = 180
COLORS = {
    "Control": "#0072B2",
    "Case": "#D55E00",
    "Accent": "#009E73",
    "Muted": "#666666",
    "Light": "#D9E4EA",
}


def _style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#d8dde2", linewidth=0.7, alpha=0.7)


def _save(fig, name):
    fig.savefig(OUT_DIR / name, dpi=DPI, facecolor="white")
    plt.close(fig)


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
    ax.set(
        xlabel="CpG methylation beta value",
        ylabel="Density",
        title="Synthetic methylation distribution",
    )
    ax.set_xlim(0, 1)
    ax.legend(frameon=False)
    _style_axis(ax)
    _save(fig, "tutorial_methylation_distribution.png")


def fragment_length_distribution(rng):
    fig, ax = plt.subplots(figsize=(6.4, 4.1), constrained_layout=True)
    for group, shift in (("Control", 0.0), ("Case", 2.0)):
        mono = rng.normal(167 + shift, 16, 48_000)
        di = rng.normal(325 + shift, 25, 5_000)
        values = np.clip(np.concatenate([mono, di]), 30, 500)
        x, density = _smooth_histogram(values, np.arange(30, 501, 2), bandwidth=3)
        ax.plot(x, 100 * density, color=COLORS[group], linewidth=2.2, label=group)
    ax.axvline(167, color="#555555", linewidth=1, linestyle="--", label="167 bp")
    ax.set(
        xlabel="Fragment length (bp)",
        ylabel="Density (scaled %)",
        title="Synthetic cfDNA fragment-length distribution",
    )
    ax.set_xlim(50, 400)
    ax.legend(frameon=False)
    _style_axis(ax)
    _save(fig, "tutorial_fragment_length_distribution.png")


def project_layout():
    fig, ax = plt.subplots(figsize=(8.2, 4.5))
    ax.axis("off")
    ax.text(
        0.03,
        0.93,
        "Example project after cftk init and cftk run",
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        color="#1f2933",
    )
    project_tree = """example_study/
|-- cftk_init.json
|-- cftk.lock.json
|-- samples.tsv
`-- results/  (created by cftk run)
    |-- 1_process/
    |-- 2_qc/
    `-- provenance/runs/"""
    ref_tree = """reference root/
`-- twist_human_methylome/
    `-- hg38/
        |-- profile.json
        |-- genome.fa
        `-- covered_targets.bed"""
    ax.text(
        0.06,
        0.78,
        project_tree,
        transform=ax.transAxes,
        family="monospace",
        fontsize=9.4,
        va="top",
        color="#1f2933",
        linespacing=1.45,
    )
    ax.text(
        0.56,
        0.78,
        ref_tree,
        transform=ax.transAxes,
        family="monospace",
        fontsize=9.4,
        va="top",
        color="#1f2933",
        linespacing=1.45,
    )
    ax.plot([0.5, 0.5], [0.13, 0.81], transform=ax.transAxes, color="#c9d5dc", lw=1)
    ax.text(
        0.06,
        0.12,
        "Project files are portable; reference components live under one selected root.",
        transform=ax.transAxes,
        fontsize=9,
        color="#4b5563",
    )
    _save(fig, "tutorial_project_layout.png")


def dinucleotide_frequency(rng):
    motifs = ["".join(pair) for pair in product("ACGT", repeat=2)]
    base = rng.dirichlet(np.full(len(motifs), 7.0))
    shifted = np.clip(base + rng.normal(0, 0.004, len(motifs)), 0.004, None)
    shifted /= shifted.sum()
    x = np.arange(len(motifs))
    fig, ax = plt.subplots(figsize=(8.1, 4.1), constrained_layout=True)
    ax.bar(x - 0.19, 100 * base, width=0.38, color=COLORS["Control"], label="Control")
    ax.bar(x + 0.19, 100 * shifted, width=0.38, color=COLORS["Case"], label="Case")
    ax.set(
        xticks=x,
        xticklabels=motifs,
        xlabel="Dinucleotide at fragment center",
        ylabel="Frequency (%)",
        title="Synthetic dinucleotide-frequency output",
    )
    ax.legend(frameon=False, ncol=2)
    _style_axis(ax)
    _save(fig, "tutorial_dinucleotide_frequency.png")


def differential_outputs(rng):
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.2))
    ax = axes[0, 0]
    control = rng.normal(loc=(-0.8, -0.1), scale=(0.45, 0.45), size=(5, 2))
    case = rng.normal(loc=(0.9, 0.25), scale=(0.45, 0.45), size=(5, 2))
    ax.scatter(control[:, 0], control[:, 1], s=40, color=COLORS["Control"], label="Control")
    ax.scatter(case[:, 0], case[:, 1], s=40, color=COLORS["Case"], label="Case")
    ax.set(xlabel="PC1 (illustrative)", ylabel="PC2 (illustrative)", title="PCA")
    ax.legend(frameon=False, fontsize=8)
    _style_axis(ax)

    ax = axes[0, 1]
    effect = rng.normal(0, 0.95, 260)
    qvalue = np.clip(rng.beta(0.8, 4.5, 260), 1e-5, 1)
    selected = (np.abs(effect) > 1.1) & (qvalue < 0.08)
    ax.scatter(effect[~selected], -np.log10(qvalue[~selected]), s=10, color="#9aa5ad", alpha=0.65)
    ax.scatter(effect[selected], -np.log10(qvalue[selected]), s=11, color=COLORS["Case"], alpha=0.8)
    ax.axhline(-np.log10(0.05), color="#777777", linestyle="--", linewidth=0.8)
    ax.set(xlabel="Mean difference", ylabel="-log10(q value)", title="DMR volcano")
    _style_axis(ax)

    ax = axes[1, 0]
    violin_values = []
    for sample_index in range(10):
        shift = 0.03 if sample_index >= 5 else 0.0
        violin_values.append(rng.beta(2.3 + shift * 8, 4.0 - shift * 4, 130))
    violins = ax.violinplot(violin_values, showmeans=False, showmedians=True)
    for index, body in enumerate(violins["bodies"]):
        body.set_facecolor(COLORS["Control"] if index < 5 else COLORS["Case"])
        body.set_edgecolor("none")
        body.set_alpha(0.72)
    ax.set(
        xticks=[1, 5, 6, 10],
        xticklabels=["C1", "C5", "K1", "K5"],
        xlabel="Illustrative samples",
        ylabel="Feature value",
        title="Feature distribution",
    )
    _style_axis(ax)

    ax = axes[1, 1]
    data = rng.normal(0, 0.55, size=(24, 10))
    data[:8, 5:] += 0.85
    data[8:16, :5] += 0.65
    image = ax.imshow(data, aspect="auto", cmap="RdBu_r", vmin=-1.8, vmax=1.8)
    ax.set(
        xticks=[2, 7],
        xticklabels=["Control", "Case"],
        yticks=[],
        title="Top-feature heatmap",
    )
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="Scaled value")
    fig.suptitle("Synthetic illustrative differential-analysis outputs", fontsize=14, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save(fig, "tutorial_differential_outputs.png")


def fragmentomics_outputs(rng):
    fig, axes = plt.subplots(2, 3, figsize=(10.6, 6.5))
    axes = axes.ravel()
    x = np.linspace(-1000, 1000, 240)

    ax = axes[0]
    occupancy = 1.0 + 0.18 * np.exp(-(x / 280) ** 2) + 0.03 * np.sin(x / 90)
    ax.plot(x, occupancy, color=COLORS["Accent"], linewidth=2)
    ax.axvline(0, color="#777777", linewidth=0.8, linestyle="--")
    ax.set(xlabel="Distance from TSS (bp)", ylabel="Normalized signal", title="Occupancy")
    _style_axis(ax)

    ax = axes[1]
    wps = 10 * np.cos(x / 145) * np.exp(-np.abs(x) / 850) + rng.normal(0, 0.6, len(x))
    ax.plot(x, wps, color=COLORS["Control"], linewidth=1.5)
    ax.axhline(0, color="#777777", linewidth=0.8)
    ax.set(xlabel="Position (bp)", ylabel="WPS", title="Window protection score")
    _style_axis(ax)

    ax = axes[2]
    bins = np.arange(1, 13)
    ratios = 0.86 + 0.13 * np.sin(bins * 0.9) + rng.normal(0, 0.025, len(bins))
    ax.bar(bins, ratios, color=COLORS["Case"], width=0.75)
    ax.axhline(1, color="#777777", linestyle="--", linewidth=0.8)
    ax.set(xlabel="Genome bin", ylabel="Short/long ratio", title="DELFI-style profile")
    _style_axis(ax)

    ax = axes[3]
    motifs = ["CCCA", "CCAG", "TGCC", "AATT", "GGCG", "TTAA"]
    values = np.sort(rng.uniform(0.04, 0.14, len(motifs)))[::-1]
    ax.bar(np.arange(len(motifs)), 100 * values, color=COLORS["Accent"])
    ax.set(xticks=np.arange(len(motifs)), xticklabels=motifs, ylabel="Frequency (%)", title="End motifs")
    _style_axis(ax)

    ax = axes[4]
    cleavage = 1.0 - 0.35 * np.exp(-(x / 135) ** 2) + 0.025 * np.cos(x / 45)
    ax.plot(x, cleavage, color=COLORS["Case"], linewidth=2)
    ax.axvline(0, color="#777777", linewidth=0.8, linestyle="--")
    ax.set(xlabel="Distance from CTCF (bp)", ylabel="Cleavage signal", title="CTCF cleavage")
    _style_axis(ax)

    axes[5].axis("off")
    fig.suptitle("Synthetic illustrative fragmentomics output types", fontsize=14, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save(fig, "tutorial_fragmentomics_outputs.png")


def mesa_outputs(rng):
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.7))
    ax = axes[0]
    fpr = np.linspace(0, 1, 120)
    ax.plot(fpr, fpr, color="#888888", linestyle="--", linewidth=0.9, label="Chance")
    ax.plot(fpr, fpr ** 0.54, color=COLORS["Accent"], linewidth=2, label="Multimodal")
    ax.plot(fpr, fpr ** 0.72, color=COLORS["Control"], linewidth=1.4, label="CpG")
    ax.set(xlabel="False-positive rate", ylabel="True-positive rate", title="LOOCV ROC")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    _style_axis(ax)

    ax = axes[1]
    correlations = np.array([[1.0, 0.46, 0.31], [0.46, 1.0, 0.38], [0.31, 0.38, 1.0]])
    image = ax.imshow(correlations, cmap="YlGnBu", vmin=0, vmax=1)
    ax.set(
        xticks=range(3),
        xticklabels=["CpG", "WPS", "DELFI"],
        yticks=range(3),
        yticklabels=["CpG", "WPS", "DELFI"],
        title="Prediction correlation",
    )
    for row in range(3):
        for column in range(3):
            ax.text(column, row, f"{correlations[row, column]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)

    ax = axes[2]
    control = np.clip(rng.normal(0.33, 0.13, 8), 0, 1)
    case = np.clip(rng.normal(0.68, 0.13, 8), 0, 1)
    ax.scatter(rng.normal(0, 0.035, len(control)), control, color=COLORS["Control"], s=28, label="Control")
    ax.scatter(rng.normal(1, 0.035, len(case)), case, color=COLORS["Case"], s=28, label="Case")
    ax.set(xticks=[0, 1], xticklabels=["Control", "Case"], ylabel="LOOCV score", title="Per-sample predictions")
    ax.legend(frameon=False, fontsize=7)
    _style_axis(ax)
    fig.suptitle("Synthetic illustrative MESA outputs", fontsize=13, y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    _save(fig, "tutorial_mesa_outputs.png")


def model_power_api(rng):
    sample_sizes = np.arange(20, 121, 10)
    auc = 0.55 + 0.29 * (1 - np.exp(-(sample_sizes - 10) / 55))
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9))
    ax = axes[0]
    ax.plot(sample_sizes, auc, marker="o", color=COLORS["Accent"], linewidth=2)
    ax.axhline(0.70, color="#777777", linewidth=0.9, linestyle="--", label="Target AUC")
    ax.set(xlabel="Total sample size", ylabel="Mean CV AUC", ylim=(0.5, 0.9), title="Power-curve summary")
    ax.legend(frameon=False, fontsize=8)
    _style_axis(ax)

    ax = axes[1]
    stability = 0.26 + 0.56 * (1 - np.exp(-(sample_sizes - 10) / 50))
    uncertainty = 0.07 * np.exp(-(sample_sizes - 20) / 90)
    ax.plot(sample_sizes, stability, marker="o", color=COLORS["Control"], linewidth=2, label="Feature recall")
    ax.fill_between(sample_sizes, stability - uncertainty, stability + uncertainty, color=COLORS["Light"], alpha=0.8)
    ax.set(xlabel="Total sample size", ylabel="Simulation diagnostic", ylim=(0, 1), title="Replicate diagnostic")
    ax.legend(frameon=False, fontsize=8)
    _style_axis(ax)
    fig.suptitle("Synthetic illustrative model-development power API output", fontsize=12, y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    _save(fig, "tutorial_model_power_api.png")


def model_power_calculator():
    sample_sizes = np.arange(20, 121, 10)
    detection = 0.20 + 0.72 * (1 - np.exp(-(sample_sizes - 12) / 48))
    attainment = 0.12 + 0.76 * (1 - np.exp(-(sample_sizes - 12) / 63))
    success = np.minimum(detection, attainment) - 0.05
    fig, ax = plt.subplots(figsize=(7.3, 4.3), constrained_layout=True)
    ax.plot(sample_sizes, detection, marker="o", color=COLORS["Control"], linewidth=2, label="Detection power")
    ax.plot(sample_sizes, attainment, marker="o", color=COLORS["Case"], linewidth=2, label="Target attainment")
    ax.plot(sample_sizes, success, marker="o", color=COLORS["Accent"], linewidth=2.3, label="Probability of success")
    ax.set(
        xlabel="Total sample size",
        ylabel="Estimated probability",
        ylim=(0, 1),
        title="Synthetic calculator-style power curve",
    )
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    _style_axis(ax)
    _save(fig, "tutorial_model_power_calculator.png")


def main():
    rng = np.random.default_rng(SEED)
    methylation_distribution(rng)
    fragment_length_distribution(rng)
    project_layout()
    dinucleotide_frequency(rng)
    differential_outputs(rng)
    fragmentomics_outputs(rng)
    mesa_outputs(rng)
    model_power_api(rng)
    model_power_calculator()


if __name__ == "__main__":
    main()
