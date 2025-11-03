import argparse
import os
import re
from typing import Any
import numpy as np
import powerlaw
import matplotlib.pyplot as plt

PI_RADIUS_SQUARED = np.pi * 0.195**2  # for φ calculation


def parse_collision_file(filename: str) -> list[float]:
    times: list[float] = []
    with open(filename, "r") as f:
        for line in f:
            if line.strip():
                t, _ = line.strip().split()
                times.append(float(t))
    if len(times) < 2:
        return []
    return [times[i + 1] - times[i] for i in range(len(times) - 1)]


def generate_semiparametric_random(fit: powerlaw.Fit, n: int) -> np.ndarray:
    data_original = np.array(fit.data_original)
    power_law = fit.power_law
    xmin = fit.xmin
    tail_data = data_original[data_original >= xmin]
    body_data = data_original[data_original < xmin]
    ntail = len(tail_data)
    n_total = len(data_original)
    p_tail = ntail / n_total

    synthetic = np.empty(n)
    for i in range(n):
        if np.random.random() < p_tail:
            synthetic[i] = power_law.generate_random(n=1)[0]
        else:
            synthetic[i] = np.random.choice(body_data)
    return synthetic


def analyze_powerlaw(
    taus: list[float], n_reps: int = 200
) -> tuple[float, float, float, float, tuple[Any, Any], tuple[Any, Any]]:
    data = np.array(taus)
    fit = powerlaw.Fit(data=data, verbose=False)
    alpha: float = fit.alpha
    alpha_std: float = fit.sigma
    xmin: float = fit.xmin or 0
    D_emp = fit.power_law.KS()

    synthetic_KS = []
    for _ in range(n_reps):
        synthetic = generate_semiparametric_random(fit, n=len(data))
        fit_syn = powerlaw.Fit(synthetic, xmin=None, verbose=False)
        D_syn = fit_syn.power_law.KS()
        synthetic_KS.append(D_syn)

    p_value = np.mean(np.array(synthetic_KS) > D_emp)

    ln_comp = fit.distribution_compare("power_law", "lognormal")
    # print(f"R = {R_ln:.3f}, p = {p_ln:.3f}")

    exp_comp = fit.distribution_compare("power_law", "exponential")

    return alpha, alpha_std, xmin, p_value, ln_comp, exp_comp


def parse_subdir_name(name: str) -> tuple[int, float]:
    """Extract N and L from a subdir name of the form 'N_%d_L%f'."""
    match = re.match(r"N_(\d+)_L([0-9.]+)", name)
    if not match:
        raise ValueError(f"Cannot parse N and L from subdir {name}")
    N = int(match.group(1))
    L = float(match.group(2))
    return N, L


def analyze_subdir(dir: str, n_reps: int):
    pattern = re.compile(r".*collisions\.csv$")
    all_taus: list[float] = []
    for fname in os.listdir(dir):
        if pattern.match(fname):
            path = os.path.join(dir, fname)
            taus = parse_collision_file(path)
            if taus:
                all_taus.extend(taus)

    if not all_taus:
        return None

    return analyze_powerlaw(all_taus, n_reps=n_reps)


def analyze_all_subdirs(parent_dir: str, n_reps: int = 200):
    """Parse subdirectories, analyze power-law for each, and return sorted arrays by φ."""

    phi_list = []
    alpha_list = []
    alpha_std_list = []
    xmin_list = []
    p_value_list = []
    ln_cmp_list = []
    exp_cmp_list = []

    for subdir in os.listdir(parent_dir):
        subdir_path = os.path.join(parent_dir, subdir)
        if not os.path.isdir(subdir_path):
            continue
        try:
            N, L = parse_subdir_name(subdir)
        except ValueError:
            continue

        analysis = analyze_subdir(subdir_path, n_reps)

        if not analysis:
            continue

        alpha, alpha_std, xmin, p_value, ln_cmp, exp_cmp = analysis
        phi = N * PI_RADIUS_SQUARED / L**2

        phi_list.append(phi)
        alpha_list.append(alpha)
        alpha_std_list.append(alpha_std)
        xmin_list.append(xmin)
        p_value_list.append(p_value)
        ln_cmp_list.append(ln_cmp)
        exp_cmp_list.append(exp_cmp)

    # Convert to NumPy arrays
    phi_arr = np.array(phi_list)
    alpha_arr = np.array(alpha_list)
    alpha_std_arr = np.array(alpha_std_list)
    xmin_arr = np.array(xmin_list)
    p_value_arr = np.array(p_value_list)
    ln_cmp_arr = np.array(ln_cmp_list)
    exp_cmp_arr = np.array(exp_cmp_list)

    # Sort by φ
    sort_idx = np.argsort(phi_arr)
    return (
        phi_arr[sort_idx],
        alpha_arr[sort_idx],
        alpha_std_arr[sort_idx],
        xmin_arr[sort_idx],
        p_value_arr[sort_idx],
        ln_cmp_arr[sort_idx],
        exp_cmp_arr[sort_idx],
    )


def plot_alpha_vs_phi(
    phi_list: np.ndarray,
    alpha_list: np.ndarray,
    alpha_std_list: np.ndarray,
    show: bool = True,
    save_path: str | None = None,
):
    """Plot α vs φ with error bars."""
    plt.errorbar(phi_list, alpha_list, yerr=alpha_std_list, fmt="o-", capsize=5)
    plt.xlabel("φ")
    plt.ylabel("α")
    plt.title("Power-law exponent α vs packing fraction φ")
    plt.grid(True)
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to {save_path}")
    if show:
        plt.show()
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Analyze power-law α vs φ for directories of simulations"
    )
    parser.add_argument(
        "directory", type=str, help="Parent directory containing subdirectories"
    )
    parser.add_argument(
        "--reps", type=int, default=200, help="Number of synthetic datasets per N"
    )
    parser.add_argument(
        "--out-dir", type=str, default=None, help="Optional path to save the plot"
    )
    args = parser.parse_args()

    if not os.path.isdir(args.directory):
        print(f"Error: {args.directory} is not a directory")
        return

    phi, alpha, alpha_std, xmin, p_value, ln_cmp, exp_cmp = analyze_all_subdirs(
        args.directory, n_reps=args.reps
    )

    for φ, a, a_std, x, p, ln, exp in zip(phi, alpha, alpha_std, xmin, p_value, ln_cmp, exp_cmp):
        print(f"φ={φ:.3f}: α={a:.3f}±{a_std:.3f}, xmin={x:.3f}, p={p:.3f}, ln_cmp=(R: {ln[0]}, p: {ln[1]}), exp_cmp=(R: {exp[0]}, p: {exp[1]})")

    plot_alpha_vs_phi(phi, alpha, alpha_std, show=True, save_path=args.out_dir)


if __name__ == "__main__":
    main()
