#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

def parse_collision_file(path):
    """Lee el archivo de colisiones con formato:
       0.0210 43
       0.0410 66
       47.7228 857 862
       Retorna lista de (t, [ids])"""
    events = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            t = float(parts[0])
            ids = list(map(int, parts[1:]))
            events.append((t, ids))
    return events

def build_unique_contact_curve(events):
    """Construye tiempos y contactos únicos acumulados."""
    seen = set()
    times, unique_counts = [], []
    total = 0
    for t, ids in sorted(events, key=lambda e: e[0]):
        new = [pid for pid in ids if pid not in seen]
        if new:
            seen.update(new)
            total += len(new)
            times.append(t)
            unique_counts.append(total)
    return np.array(times), np.array(unique_counts)


def compute_slope(time, contacts, t_start=10.0, min_points=5):
    """Ajusta una recta a partir de t >= t_start (por defecto 10 s)."""
    time = np.asarray(time)
    contacts = np.asarray(contacts)

    if len(time) < min_points:
        return np.nan, np.nan

    # primer índice con t >= t_start
    n0 = np.searchsorted(time, t_start, side='left')
    
    if len(time) - n0 < min_points:
        return np.nan, np.nan

    slope, _, _, _, std_err = linregress(time[n0:], contacts[n0:])
    return slope, std_err

def extract_phi(dirname):
    """Calcula ϕ ≈ N·πr²/L² a partir del nombre del directorio N_XXX_LYYY."""
    try:
        parts = dirname.split("_")
        N = int(parts[1])
        L = float(parts[2][1:])
        r_avg = 0.195
        return N * np.pi * (r_avg ** 2) / (L ** 2)
    except Exception:
        return np.nan

def process_all(base_dir="outputs"):
    results = []
    for subdir in sorted(os.listdir(base_dir)):
        folder = os.path.join(base_dir, subdir)
        if not os.path.isdir(folder):
            continue

        files = [f for f in os.listdir(folder) if "collision" in f or "collisions" in f]
        if not files:
            print(f"⚠️ No se encontró archivo de colisiones en {subdir}")
            continue

        slopes = []
        for f in files:
            try:
                events = parse_collision_file(os.path.join(folder, f))
                time, contacts = build_unique_contact_curve(events)
                slope, _ = compute_slope(time, contacts, t_start=10.0)
                slopes.append(slope)
                plt.plot(time, contacts, label=f"{subdir}")
            except Exception as e:
                print(f"Error procesando {f}: {e}")

        if slopes:
            phi = extract_phi(subdir)
            Q = np.mean(slopes)
            Qerr = np.std(slopes, ddof=0)  # <- desvío estándar absoluto
            results.append((phi, Q, Qerr))
            print(f"{subdir}: ϕ={phi:.3f}, Q={Q:.4f} ± {Qerr:.4f}")

    plt.xlabel("Tiempo [s]")
    plt.ylabel("Contactos únicos acumulados")
    plt.title("Curvas de contactos únicos (SFM)")
    plt.legend(fontsize=8)
    plt.tight_layout()
    os.makedirs("images", exist_ok=True)
    # plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.subplots_adjust(top=0.93, bottom=0.12, left=0.12, right=0.95)
    plt.savefig("images/curves_contacts.png", dpi=300)

    # plt.savefig("images/curves_contacts.png", dpi=300)
    plt.close()

    return sorted(results, key=lambda x: x[0])

def plot_Q_vs_phi(results, output="images/Q_vs_phi.png"):
    phi, Q, Qerr = np.array(results).T
    plt.errorbar(phi, Q, yerr=Qerr, fmt='o-', capsize=4)
    plt.xlabel(r"Fracción de área ocupada $\phi$", fontsize=14)
    plt.ylabel(r"Pendiente $\mathit{Q}$ (s$^{-1}$)", fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    print(f"📈 Figura guardada en {output}")

if __name__ == "__main__":
    results = process_all("PedestrianDynamic/OutputsOriginales")
    if results:
        plot_Q_vs_phi(results)
