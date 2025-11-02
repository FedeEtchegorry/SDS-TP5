import argparse
import os
import re
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np

import powerlaw


def parse_collision_file(filename: str) -> List[float]:
    """Reads a file with lines '<t> <id>' and returns list of time intervals tau_i."""
    times: List[float] = []
    with open(filename, "r") as f:
        for line in f:
            if line.strip():
                t, _ = line.strip().split()
                times.append(float(t))

    if len(times) < 2:
        return []  # skip files with not enough data

    taus: List[float] = [times[i + 1] - times[i] for i in range(len(times) - 1)]
    return taus


def plot_collision_intervals(taus: List[float], out_file: str) -> None:
    """Plots tau_i and saves the figure to out_file."""
    plt.figure(figsize=(8, 4))
    plt.plot(taus, marker="o", linestyle="-", color="blue")
    plt.xlabel("Collision index")
    plt.ylabel(r"Time between collisions $\tau_i$")
    plt.title("Time Intervals Between Collisions")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    plt.close()
    print(f"Plot saved to {out_file}")


def analyze_powerlaw(taus: List[float]) -> Tuple[float, float, float, float]:
    """Analyzes tau_i for power-law distribution. Returns alpha, xmin, R, p."""
    fit = powerlaw.Fit(taus, xmin=min(taus))
    alpha: float = fit.power_law.alpha
    xmin: float = fit.power_law.xmin
    R, p = fit.distribution_compare("power_law", "lognormal")

    print(f"Estimated alpha: {alpha:.3f}")
    print(f"Estimated xmin: {xmin:.3f}")
    print(f"Comparison with lognormal: R = {R:.3f}, p = {p:.3f}")
    if R > 0:
        print("Power law is favored over lognormal")
    else:
        print("Lognormal is favored over power law")

    return alpha, xmin, R, p


def main():
    parser = argparse.ArgumentParser(
        description="Parse all collision files in a directory, merge them, and analyze power-law distribution."
    )
    parser.add_argument(
        "directory", type=str, help="Directory containing collision files"
    )
    parser.add_argument(
        "--out_dir", type=str, default="powerlaw", help="Folder to save plots"
    )
    args = parser.parse_args()

    if not os.path.isdir(args.directory):
        print(f"Error: {args.directory} is not a directory")
        return

    os.makedirs(args.out_dir, exist_ok=True)
    pattern = re.compile(r".*collisions\.csv$")
    all_taus: List[float] = []

    # Scan directory for matching files
    for fname in os.listdir(args.directory):
        if pattern.match(fname):
            path = os.path.join(args.directory, fname)
            taus = parse_collision_file(path)
            if taus:
                all_taus.extend(taus)

    if not all_taus:
        print("No valid collision data found in directory.")
        return

    # Save plot
    # out_file = os.path.join(args.out_dir, "merged_collision_intervals.png")
    # plot_collision_intervals(all_taus, out_file)

    # Analyze power-law
    # analyze_powerlaw(all_taus)
    fit = powerlaw.Fit(all_taus)
    alpha = fit.power_law.alpha
    xmin = fit.power_law.xmin
    D = fit.power_law.KS()

    # Bootstrap for uncertainty
    alphas = []
    for _ in range(1000):
        resample = np.random.choice(all_taus, len(all_taus), replace=True)
        fit_resample = powerlaw.Fit(resample, xmin=xmin)
        alphas.append(fit_resample.power_law.alpha)

    alpha_err = np.std(alphas)

    R, p = fit.distribution_compare("power_law", "lognormal")

    print(f"α = {alpha:.3f} ± {alpha_err:.3f},  xmin = {xmin:.3f},  KS = {D:.3f}, (R, p) = ({R}, {p})")


if __name__ == "__main__":
    main()
