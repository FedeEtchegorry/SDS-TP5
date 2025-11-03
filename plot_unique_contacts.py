import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from glob import glob

# ==========================================================
#  LECTURA DE ARCHIVO DE COLISIONES
# ==========================================================
def load_curve(file_path):
    """Lee un archivo de colisiones, permitiendo múltiples IDs por mismo tiempo."""
    times, counts = [], []
    unique_ids = set()

    try:
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split()
                try:
                    time = float(parts[0])
                except ValueError:
                    continue

                ids = []
                for p in parts[1:]:
                    try:
                        ids.append(int(p))
                    except ValueError:
                        pass

                unique_ids.update(ids)
                times.append(time)
                counts.append(len(unique_ids))

    except Exception as e:
        print(f"⚠️ Error leyendo {file_path}: {e}")
        return None, None

    if not times:
        print(f"⚠️ Archivo vacío o ilegible: {file_path}")
        return None, None

    return np.array(times), np.array(counts)

# ==========================================================
#  PROMEDIO DE TODAS LAS CURVAS EN UNA CARPETA
# ==========================================================
def plot_average_contacts(folder_path, output_path=None):
    """Promedia todas las curvas de una carpeta y grafica con barras de error."""
    files = sorted(glob(os.path.join(folder_path, "*collisions*.csv")))
    if not files:
        print(f"⚠️ No se encontraron archivos de colisiones en: {folder_path}")
        return

    print(f"📂 Archivos encontrados: {len(files)}")

    # Cargar todas las curvas y buscar rango temporal global
    curves = []
    t_min, t_max = float("inf"), 0.0

    for f in files:
        t, c = load_curve(f)
        if t is None or len(t) < 2:
            continue
        curves.append((t, c))
        t_min = min(t_min, t[0])
        t_max = max(t_max, t[-1])

    if not curves:
        print("⚠️ No se pudieron leer curvas válidas.")
        return

    # Eje de tiempo común (uniforme)
    num_points = max(len(t) for t, _ in curves)
    common_time = np.linspace(t_min, t_max, num_points)

    # Interpolar todas las curvas al eje común
    all_curves = []
    for t, c in curves:
        interp_c = np.interp(common_time, t, c)
        all_curves.append(interp_c)

    all_curves = np.array(all_curves)
    mean_curve = np.mean(all_curves, axis=0)
    std_curve = np.std(all_curves, axis=0)

    # ==========================================================
    #  GRAFICAR RESULTADO
    # ==========================================================
    plt.figure(figsize=(8,5))
    plt.plot(common_time, mean_curve, color='tab:blue', lw=2, label='Datos')
    plt.fill_between(common_time,
                     mean_curve - std_curve,
                     mean_curve + std_curve,
                     color='tab:blue', alpha=0.25, label='Desvío estándar')

    plt.xlabel("Tiempo (s)", fontsize=14)
    plt.ylabel(r"Contactos únicos acumulados $N_c$", fontsize=14)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.axvline(x=10, color='red', linestyle='--', linewidth=1.6)
    plt.legend(fontsize=12)
    plt.tick_params(axis='both', labelsize=12)
    plt.tight_layout()

    # Guardar o mostrar
    if output_path:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        plt.savefig(output_path, dpi=300)
        print(f"✅ Gráfico guardado en: {output_path}")
    else:
        plt.show()

# ==========================================================
#  MAIN
# ==========================================================
if __name__ == "__main__":
    if len(sys.argv) not in [2, 3]:
        print("Uso: python plot_unique_contacts.py <carpeta_N_xx_Lx.xxx> [nombre_salida.png]")
        sys.exit(1)

    folder_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) == 3 else None
    plot_average_contacts(folder_path, output_path)
