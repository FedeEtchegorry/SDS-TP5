import argparse
import os
import re
from typing import Any, List, Tuple, Optional
import numpy as np
import powerlaw
import matplotlib.pyplot as plt
from scipy.special import erfc

PI_RADIUS_SQUARED = np.pi * 0.195**2  # for φ calculation


def parse_collision_file(filename: str) -> List[float]:
    times: List[float] = []
    with open(filename, "r") as f:
        for lineno, line in enumerate(f, start=1):
            s = line.strip()
            if not s:
                continue  # línea vacía
            parts = s.split()          # puede haber 2, 3, ... tokens
            try:
                t = float(parts[0])    # SOLO el primer token es el tiempo
            except (ValueError, IndexError):
                # Si hay encabezados/comentarios, los saltamos
                # print(f"Saltando línea {lineno}: {s}")
                continue
            times.append(t)

    if len(times) < 2:
        return []
    # Δτ_i = t_{i+1} - t_i (solo positivos)
    taus = [times[i+1] - times[i] for i in range(len(times)-1)]
    # Filtramos cualquier Δτ no positivo (si hubiera duplicados)
    return [dt for dt in taus if dt > 0.0]

def generate_semiparametric_random(fit: powerlaw.Fit, n: int) -> np.ndarray:
    data_original = np.array(fit.data_original)
    power_law = fit.power_law
    xmin = fit.xmin
    tail_data = data_original[data_original >= xmin]
    body_data = data_original[data_original < xmin]
    ntail = len(tail_data)
    n_total = len(data_original)
    p_tail = ntail / n_total if n_total > 0 else 0.0

    synthetic = np.empty(n)
    for i in range(n):
        if np.random.random() < p_tail and ntail > 0:
            synthetic[i] = power_law.generate_random(n=1)[0]
        else:
            synthetic[i] = np.random.choice(body_data if len(body_data) > 0 else data_original)
    return synthetic


def analyze_powerlaw(
        taus: List[float], n_reps: int = 200
) -> Tuple[float, float, float, float, Tuple[Any, Any], Tuple[Any, Any]]:
    data = np.array(taus, dtype=float)
    data = data[data > 0]
    if data.size == 0:
        raise ValueError("No positive data for analysis.")

    fit = powerlaw.Fit(data=data, verbose=False)
    alpha: float = fit.alpha
    alpha_std: float = fit.sigma
    xmin: float = fit.xmin or 0.0
    D_emp = fit.power_law.KS()

    synthetic_KS = []
    for _ in range(n_reps):
        synthetic = generate_semiparametric_random(fit, n=len(data))
        fit_syn = powerlaw.Fit(synthetic, xmin=None, verbose=False)
        D_syn = fit_syn.power_law.KS()
        synthetic_KS.append(D_syn)

    p_value = float(np.mean(np.array(synthetic_KS) > D_emp))

    ln_comp = fit.distribution_compare("power_law", "lognormal")     # (R, p)
    exp_comp = fit.distribution_compare("power_law", "exponential")  # (R, p)

    return alpha, alpha_std, xmin, p_value, ln_comp, exp_comp


def parse_subdir_name(name: str) -> Tuple[int, float]:
    """Extract N and L from a subdir name of the form 'N_%d_L%f'."""
    match = re.match(r"N_(\d+)_L([0-9.]+)", name)
    if not match:
        raise ValueError(f"Cannot parse N and L from subdir {name}")
    N = int(match.group(1))
    L = float(match.group(2))
    return N, L


def analyze_subdir(dir_path: str, n_reps: int) -> Optional[Tuple[
    Tuple[float, float, float, float, Tuple[Any, Any], Tuple[Any, Any]],
    List[float]
]]:
    pattern = re.compile(r".*collisions\.csv$")
    all_taus: List[float] = []
    for fname in os.listdir(dir_path):
        if pattern.match(fname):
            path = os.path.join(dir_path, fname)
            taus = parse_collision_file(path)
            if taus:
                all_taus.extend(taus)

    if not all_taus:
        return None

    analysis = analyze_powerlaw(all_taus, n_reps=n_reps)
    return analysis, all_taus


def log_hist_pdf(data: np.ndarray, nbins: int = 30, ax: Optional[plt.Axes] = None):
    """PDF con bins logarítmicos (normalizado por ancho de bin)."""
    x = np.asarray(data, dtype=float)
    x = x[x > 0]
    if x.size == 0:
        return ax

    xmin, xmax = x.min(), x.max()
    if xmin <= 0 or xmin == xmax:
        return ax

    edges = np.logspace(np.log10(xmin), np.log10(xmax), nbins + 1)
    counts, edges = np.histogram(x, bins=edges)
    widths = np.diff(edges)
    pdf = counts / (counts.sum() * widths)
    centers = np.sqrt(edges[:-1] * edges[1:])

    if ax is None:
        ax = plt.gca()
    ax.loglog(centers, pdf, "o", label="PDF (bins log)")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$p(\tau)$")
    ax.grid(True, which="both", ls=":")
    return ax


def ccdf_loglog(data: np.ndarray, ax: Optional[plt.Axes] = None):
    """CCDF (1-CDF) en escala log–log."""
    x = np.sort(np.asarray(data, dtype=float))
    x = x[x > 0]
    if x.size == 0:
        return ax
    n = x.size
    ccdf = 1.0 - (np.arange(1, n + 1) / n)


    if ax is None:
        ax = plt.gca()
    ax.loglog(x, ccdf, "o", label="CCDF")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel("P(T ≥ τ)")
    ax.grid(True, which="both", ls=":")
    return ax


def semilog_pdf(data: np.ndarray, nbins: int = 30, ax: Optional[plt.Axes] = None):
    """PDF con bins lineales y eje y logarítmico (útil para exponenciales)."""
    x = np.asarray(data, dtype=float)
    x = x[x > 0]
    if x.size == 0:
        return ax
    xmin, xmax = x.min(), x.max()
    if xmin == xmax:
        return ax
    edges = np.linspace(xmin, xmax, nbins + 1)
    counts, edges = np.histogram(x, bins=edges, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])

    if ax is None:
        ax = plt.gca()
    ax.semilogy(centers, counts, "o", label="PDF (bins lineales)")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$p(\tau)$")
    ax.grid(True, which="both", ls=":")
    return ax

def plot_three_views(
        data: List[float],
        title_prefix: str,
        nbins: int = 30,
        show: bool = False,
        save_prefix: Optional[str] = None,
):
    """Genera tres figuras: PDF (log bins), CCDF (log–log) y PDF (semi-log),
    superponiendo modelos power-law (cola), exponencial y lognormal **condicionados a x≥xmin**."""
    arr = np.asarray(data, dtype=float)
    arr = arr[arr > 0]
    if arr.size == 0:
        return

    # ---------- Ajuste power-law en cola ----------
    try:
        fit_tmp = powerlaw.Fit(arr, verbose=False)
        alpha_hat = float(fit_tmp.alpha)
        xmin_hat  = float(fit_tmp.xmin)
    except Exception:
        alpha_hat = None
        xmin_hat  = None

    if (alpha_hat is None) or (xmin_hat is None) or (xmin_hat <= 0):
        return

    # ---------- Datos de cola ----------
    x_tail = arr[arr >= xmin_hat]
    if x_tail.size == 0:
        return

    # ---------- Parámetros de modelos condicionados (cola) ----------
    # Exponencial cond.: f(x|x≥xmin) = λ e^{-λ (x-xmin)},  λ = 1/(E[X|tail]-xmin)
    lam_tail = None
    m_tail = float(x_tail.mean())
    if m_tail > xmin_hat:
        lam_tail = 1.0 / (m_tail - xmin_hat)

    # Lognormal cond.: f(x|x≥xmin) = f(x)/S(xmin) con μ,σ estimados en cola
    logx_tail  = np.log(x_tail)
    mu_tail    = float(logx_tail.mean())
    sigma_tail = float(logx_tail.std(ddof=0))
    S_logn_xmin = None
    if sigma_tail > 0:
        zmin = (np.log(xmin_hat) - mu_tail) / sigma_tail
        S_logn_xmin = 0.5 * erfc(zmin / np.sqrt(2.0))  # = 1 - Φ(zmin)

    # ---------- Ejes x ----------
    x_min  = float(arr.min())
    x_max  = float(arr.max())
    if not (np.isfinite(x_min) and np.isfinite(x_max)) or x_min <= 0 or x_min == x_max:
        return

    x_plot = np.logspace(np.log10(x_min),  np.log10(x_max),  400)
    x_plot_tail = x_plot[x_plot >= xmin_hat]
    if x_plot_tail.size == 0:
        return

    # ---------- Curvas teóricas en cola (PDF y CCDF) ----------
    # Power-law (pdf y ccdf)
    mask_pl = x_plot >= xmin_hat
    pdf_pl  = np.zeros_like(x_plot)
    ccdf_pl = np.zeros_like(x_plot)
    if alpha_hat > 1:
        pdf_pl[mask_pl]  = (alpha_hat - 1.0) / xmin_hat * (x_plot[mask_pl] / xmin_hat) ** (-alpha_hat)
        ccdf_pl[mask_pl] = (x_plot[mask_pl] / xmin_hat) ** (-(alpha_hat - 1.0))

    # Exponencial condicionada
    pdf_exp_tail  = None
    ccdf_exp_tail = None
    if lam_tail is not None:
        dx_tail = x_plot_tail - xmin_hat
        pdf_exp_tail  = lam_tail * np.exp(-lam_tail * dx_tail)
        ccdf_exp_tail = np.exp(-lam_tail * dx_tail)

    # Lognormal condicionada
    pdf_logn_tail  = None
    ccdf_logn_tail = None
    if (sigma_tail > 0) and (S_logn_xmin is not None) and (S_logn_xmin > 0):
        z_tail = (np.log(x_plot_tail) - mu_tail) / sigma_tail
        base_pdf_logn  = (1.0 / (x_plot_tail * sigma_tail * np.sqrt(2.0*np.pi))) \
                         * np.exp(-0.5 * z_tail**2)
        base_ccdf_logn = 0.5 * erfc(z_tail / np.sqrt(2.0))
        pdf_logn_tail  = base_pdf_logn  / S_logn_xmin
        ccdf_logn_tail = base_ccdf_logn / S_logn_xmin  # = S(x)/S(xmin)

    # ---------- 1) PDF con bins logarítmicos + modelos (cola) ----------
    fig1, ax1 = plt.subplots()
    # PDF empírica (sobre todo el rango; es OK para ver el cuerpo)
    edges = np.logspace(np.log10(x_min), np.log10(x_max), nbins + 1)
    counts, _ = np.histogram(arr, bins=edges)
    widths     = np.diff(edges)
    centers    = np.sqrt(edges[:-1] * edges[1:])
    pdf_emp    = counts / (counts.sum() * widths)
    ax1.loglog(centers, pdf_emp, "o", label="PDF (bins log)")

    # Modelos SOLO en cola
    if pdf_exp_tail is not None:
        ax1.loglog(x_plot_tail, pdf_exp_tail, "-", label="Exp (pdf) cond. x≥xmin")
    if pdf_logn_tail is not None:
        ax1.loglog(x_plot_tail, pdf_logn_tail, "-", label="Lognormal (pdf) cond. x≥xmin")
    if alpha_hat > 1:
        ax1.loglog(x_plot[mask_pl], pdf_pl[mask_pl], "-", label=f"Power law (pdf), α={alpha_hat:.2f}, xmin={xmin_hat:.2g}")

    ax1.set_xlabel(r"$\tau$", fontsize=12)
    ax1.set_ylabel(r"$p(\tau)$", fontsize=12)
    ax1.set_title(f"{title_prefix} — PDF (bins log)")
    ax1.grid(True, which="both", ls=":")
    ax1.legend()
    if save_prefix:
        fig1.savefig(f"{save_prefix}_pdf_logbins.png", dpi=300, bbox_inches="tight")
        print(f"Saved {save_prefix}_pdf_logbins.png")
    if show:
        plt.show()
    plt.close(fig1)

    # ---------- 2) CCDF log–log SOLO COLA (condicionada a x>=xmin) ----------
    fig2, ax2 = plt.subplots()

    # CCDF empírica condicionada: P(T >= x | T >= xmin)
    xs = np.sort(arr)
    mask_tail = xs >= xmin_hat
    xs_tail = xs[mask_tail]
    n_tail = xs_tail.size
    if n_tail == 0:
        return

    ccdf_tail = 1.0 - (np.arange(1, n_tail + 1) / n_tail)  # = S(x)/S(xmin)
    ax2.loglog(xs_tail, ccdf_tail, "o", label="CCDF empírica (cola)")

    # Modelos teóricos (ya calculados como condicionados a x>=xmin)
    if ccdf_exp_tail is not None:
        ax2.loglog(x_plot_tail, ccdf_exp_tail, "-", label="Exp (ccdf) cond.")
    if ccdf_logn_tail is not None:
        ax2.loglog(x_plot_tail, ccdf_logn_tail, "-", label="Lognormal (ccdf) cond.")
    if alpha_hat > 1:
        ax2.loglog(x_plot[mask_pl], ccdf_pl[mask_pl], "-", label="Power law (ccdf)")

    ax2.set_xlabel(r"$\tau$")
    ax2.set_ylabel(r"$P(T \geq \tau)$")

    # Limitar ejes para ver solo la cola con buena lectura
    ax2.set_xlim(left=xmin_hat, right=x_max * 1.05)
    ymin = max(ccdf_tail.min() * 0.7, 1e-4)
    ax2.set_ylim(bottom=ymin, top=1.1)

    ax2.grid(True, which="both", ls=":")
    ax2.legend()
    if save_prefix:
        fig2.savefig(f"{save_prefix}_ccdf_loglog.png", dpi=300, bbox_inches="tight")
        print(f"Saved {save_prefix}_ccdf_loglog.png")
    if show:
        plt.show()
    plt.close(fig2)

    # ---------- 3) PDF semi-log + modelos (cola) ----------
    fig3, ax3 = plt.subplots()
    edges_lin   = np.linspace(x_min, x_max, nbins + 1)
    counts_lin, edges_lin = np.histogram(arr, bins=edges_lin, density=True)
    centers_lin = 0.5 * (edges_lin[:-1] + edges_lin[1:])
    ax3.semilogy(centers_lin, counts_lin, "o", label="PDF (bins lineales)")

    if pdf_exp_tail is not None:
        ax3.semilogy(x_plot_tail, pdf_exp_tail, "-", label="Exp (pdf) cond.")
    if pdf_logn_tail is not None:
        ax3.semilogy(x_plot_tail, pdf_logn_tail, "-", label="Lognormal (pdf) cond.")
    if alpha_hat > 1:
        ax3.semilogy(x_plot[mask_pl], pdf_pl[mask_pl], "-", label="Power law (pdf)")

    ax3.set_xlabel(r"$\tau$", fontsize=12)
    ax3.set_ylabel(r"$p(\tau)$", fontsize=12)

    ax3.grid(True, which="both", ls=":")
    ax3.legend()
    if save_prefix:
        fig3.savefig(f"{save_prefix}_pdf_semilog.png", dpi=300, bbox_inches="tight")
        print(f"Saved {save_prefix}_pdf_semilog.png")
    if show:
        plt.show()
    plt.close(fig3)


def analyze_all_subdirs(parent_dir: str, n_reps: int = 200):
    """
    Parse subdirectories, analyze power-law for each, and return sorted arrays by φ,
    junto con listas auxiliares para plotting por subdir.
    """
    phi_list = []
    alpha_list = []
    alpha_std_list = []
    xmin_list = []
    p_value_list = []
    ln_cmp_list = []
    exp_cmp_list = []

    subdir_names: List[str] = []
    subdir_taus: List[List[float]] = []

    for subdir in os.listdir(parent_dir):
        subdir_path = os.path.join(parent_dir, subdir)
        if not os.path.isdir(subdir_path):
            continue
        try:
            N, L = parse_subdir_name(subdir)
        except ValueError:
            continue

        result = analyze_subdir(subdir_path, n_reps)
        if not result:
            continue

        analysis, taus = result
        alpha, alpha_std, xmin, p_value, ln_cmp, exp_cmp = analysis
        phi = N * PI_RADIUS_SQUARED / L**2

        phi_list.append(phi)
        alpha_list.append(alpha)
        alpha_std_list.append(alpha_std)
        xmin_list.append(xmin)
        p_value_list.append(p_value)
        ln_cmp_list.append(ln_cmp)
        exp_cmp_list.append(exp_cmp)

        subdir_names.append(subdir)
        subdir_taus.append(taus)

    # Convert to NumPy arrays
    phi_arr = np.array(phi_list)
    alpha_arr = np.array(alpha_list)
    alpha_std_arr = np.array(alpha_std_list)
    xmin_arr = np.array(xmin_list)
    p_value_arr = np.array(p_value_list)
    ln_cmp_arr = np.array(ln_cmp_list, dtype=object)
    exp_cmp_arr = np.array(exp_cmp_list, dtype=object)

    # Sort by φ (and permute acompañantes)
    sort_idx = np.argsort(phi_arr)
    phi_arr = phi_arr[sort_idx]
    alpha_arr = alpha_arr[sort_idx]
    alpha_std_arr = alpha_std_arr[sort_idx]
    xmin_arr = xmin_arr[sort_idx]
    p_value_arr = p_value_arr[sort_idx]
    ln_cmp_arr = ln_cmp_arr[sort_idx]
    exp_cmp_arr = exp_cmp_arr[sort_idx]
    subdir_names = [subdir_names[i] for i in sort_idx]
    subdir_taus = [subdir_taus[i] for i in sort_idx]

    return (
        phi_arr,
        alpha_arr,
        alpha_std_arr,
        xmin_arr,
        p_value_arr,
        ln_cmp_arr,
        exp_cmp_arr,
        subdir_names,
        subdir_taus,
    )


def plot_alpha_vs_phi(
        phi_list: np.ndarray,
        alpha_list: np.ndarray,
        alpha_std_list: np.ndarray,
        show: bool = True,
        save_path: Optional[str] = None,
):
    """Plot α vs φ with error bars."""
    plt.errorbar(phi_list, alpha_list, yerr=alpha_std_list, fmt="o-", capsize=5)
    plt.xlabel(r"$\varphi$")
    plt.ylabel(r"$\alpha$")
    plt.grid(True)
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to {save_path}")
    if show:
        plt.show()
    plt.close()


def sanitize(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze power-law α vs φ for directories of simulations, with optional histograms/CCDF."
    )
    parser.add_argument("directory", type=str, help="Parent directory containing subdirectories")
    parser.add_argument("--reps", type=int, default=200, help="Number of synthetic datasets per N")
    parser.add_argument("--out-dir", type=str, default=None, help="Optional path to save plots")
    parser.add_argument("--plots", action="store_true", help="Generate PDF/CCDF plots per subdir")
    parser.add_argument("--nbins", type=int, default=30, help="Number of bins for histograms")
    parser.add_argument("--no-show", action="store_true", help="Do not display plots interactively")

    parser.add_argument("--overlay-datasets", action="store_true",
                    help="Superponer todas las corridas en una sola figura (CCDF y PDF log-bins).")
    args = parser.parse_args()
    if not os.path.isdir(args.directory):
        print(f"Error: {args.directory} is not a directory")
        return

    (
        phi, alpha, alpha_std, xmin, p_value, ln_cmp, exp_cmp, subdir_names, subdir_taus
    ) = analyze_all_subdirs(args.directory, n_reps=args.reps)

    for φ, a, a_std, x, p, ln, exp, name in zip(
            phi, alpha, alpha_std, xmin, p_value, ln_cmp, exp_cmp, subdir_names
    ):
        print(
            f"φ={φ:.3f}: α={a:.3f}±{a_std:.3f}, xmin={x:.3f}, p={p:.3f}, "
            f"ln_cmp=(R: {ln[0]:.4g}, p: {ln[1]:.3g}), "
            f"exp_cmp=(R: {exp[0]:.4g}, p: {exp[1]:.3g})  [{name}]"
        )

    # Plot α vs φ
    out_alpha_phi = None
    if args.out_dir:
        os.makedirs(args.out_dir, exist_ok=True)
        out_alpha_phi = os.path.join(args.out_dir, "alpha_vs_phi.png")

    plot_alpha_vs_phi(phi, alpha, alpha_std, show=not args.no_show, save_path=out_alpha_phi)

    # Optionally produce three-view plots per subdir
    if args.plots:
        for name, taus, φ in zip(subdir_names, subdir_taus, phi):
            title = f"{name} (φ={φ:.3f})"
            save_prefix = None
            if args.out_dir:
                base = f"{sanitize(name)}__phi_{φ:.3f}"
                save_prefix = os.path.join(args.out_dir, base)
            plot_three_views(
                taus,
                title_prefix=title,
                nbins=args.nbins,
                show=not args.no_show,
                save_prefix=save_prefix,
            )
    if args.overlay_datasets and subdir_taus:
        # Un color por dataset (matplotlib cicla colores)
        # 2.a) CCDF log-log superpuesta
        fig_ccdf, ax_ccdf = plt.subplots()
        for name, taus, φ in zip(subdir_names, subdir_taus, phi):
            arr = np.asarray(taus, dtype=float)
            arr = arr[arr > 0]
            if arr.size == 0:
                continue
            xs = np.sort(arr)
            n = xs.size
            ccdf = 1.0 - (np.arange(1, n + 1) / n)
            ax_ccdf.loglog(xs, ccdf, '.', ms=3, alpha=0.7, label=f"{name} (φ={φ:.3f})")
        ax_ccdf.set_xlabel(r"$\tau$")
        ax_ccdf.set_ylabel("P(T ≥ τ)")
        ax_ccdf.set_title("CCDF (todas las corridas)")
        ax_ccdf.grid(True, which="both", ls=":")
        ax_ccdf.legend(fontsize=8, ncol=2)
        if args.out_dir:
            out_path = os.path.join(args.out_dir, "ALL_ccdf_loglog.png")
            fig_ccdf.savefig(out_path, dpi=300, bbox_inches="tight")
            print(f"Saved {out_path}")
        if not args.no_show:
            plt.show()
        plt.close(fig_ccdf)

        # 2.b) PDF con bins logarítmicos superpuesta (mismos bins globales)
        # Definimos bins comunes usando el min/max global:
        all_vals = np.concatenate([np.asarray(t, float) for t in subdir_taus])
        all_vals = all_vals[all_vals > 0]
        if all_vals.size > 0:
            xmin_g, xmax_g = float(all_vals.min()), float(all_vals.max())
            if xmin_g > 0 and xmin_g < xmax_g:
                edges_g = np.logspace(np.log10(xmin_g), np.log10(xmax_g), args.nbins + 1)

                fig_pdf, ax_pdf = plt.subplots()
                for name, taus, φ in zip(subdir_names, subdir_taus, phi):
                    arr = np.asarray(taus, dtype=float)
                    arr = arr[arr > 0]
                    if arr.size == 0:
                        continue
                    counts, _ = np.histogram(arr, bins=edges_g)
                    widths = np.diff(edges_g)
                    pdf = counts / (counts.sum() * widths) if counts.sum() > 0 else counts.astype(float)
                    centers = np.sqrt(edges_g[:-1] * edges_g[1:])
                    ax_pdf.loglog(centers, pdf, '-', lw=1.2, alpha=0.8, label=f"{name} (φ={φ:.3f})")

                ax_pdf.set_xlabel(r"$\tau$")
                ax_pdf.set_ylabel(r"$p(\tau)$")
                ax_pdf.grid(True, which="both", ls=":")
                ax_pdf.legend(fontsize=8, ncol=2)
                if args.out_dir:
                    out_path = os.path.join(args.out_dir, "ALL_pdf_logbins.png")
                    fig_pdf.savefig(out_path, dpi=300, bbox_inches="tight")
                    print(f"Saved {out_path}")
                if not args.no_show:
                    plt.show()
                plt.close(fig_pdf)

if __name__ == "__main__":
    main()
