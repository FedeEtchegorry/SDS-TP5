import argparse
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
import numpy as np


def read_simulation_output(path):
    steps = []
    with open(path, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        time = float(lines[i].strip())
        i += 1
        particles = []
        while i < len(lines) and ',' in lines[i]:
            parts = lines[i].strip().split(',')
            if len(parts) < 6:
                break
            pid = int(parts[0])
            x, y, vx, vy, r = map(float, parts[1:])
            particles.append((pid, x, y, vx, vy, r))
            i += 1
        steps.append((time, particles))
    return steps


def animate_simulation(file_path, L=6.0, save=False, output_name="animation.mp4"):
    steps = read_simulation_output(file_path)
    fig, ax = plt.subplots()
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_aspect('equal')
    ax.set_title("Pedestrian Simulation (SFM)")

    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

    circles = {}
    labels = {}
    collided = set()  # IDs que ya colisionaron alguna vez

    # Inicialización
    for pid, x, y, vx, vy, r in steps[0][1]:
        color = 'red' if pid == 0 else 'royalblue'
        c = Circle((x, y), r, color=color, alpha=0.7)
        ax.add_patch(c)
        circles[pid] = c
        label = ax.text(x, y, str(pid),
                        ha='center', va='center',
                        fontsize=6, color='white', weight='bold')
        labels[pid] = label

    def init():
        time_text.set_text('')
        return list(circles.values()) + list(labels.values()) + [time_text]

    def update(frame):
        time, particles = steps[frame]
        current_ids = {p[0] for p in particles}
        x0, y0, r0 = particles[0][1], particles[0][2], particles[0][5]

        # 🔹 borrar fantasmas (IDs que ya no existen)
        for old_id in list(circles.keys()):
            if old_id not in current_ids:
                circles[old_id].remove()
                labels[old_id].remove()
                del circles[old_id]
                del labels[old_id]
                collided.discard(old_id)

        # 🔹 actualizar o crear partículas
        for pid, x, y, vx, vy, r in particles:
            if pid not in circles:
                color = 'red' if pid == 0 else 'royalblue'
                c = Circle((x, y), r, color=color, alpha=0.7)
                ax.add_patch(c)
                circles[pid] = c
                label = ax.text(x, y, str(pid),
                                ha='center', va='center',
                                fontsize=6, color='white', weight='bold')
                labels[pid] = label

            c = circles[pid]
            label = labels[pid]
            c.center = (x, y)
            label.set_position((x, y))

            # contacto con centro → marca como colisionada
            if pid != 0:
                dist = np.hypot(x - x0, y - y0)
                if dist < (r + r0):
                    collided.add(pid)

            # actualizar color
            if pid in collided or pid == 0:
                c.set_color('red')
            else:
                c.set_color('royalblue')

        time_text.set_text(f"t = {time:.2f} s")
        return list(circles.values()) + list(labels.values()) + [time_text]

    ani = animation.FuncAnimation(
        fig, update, frames=len(steps),
        init_func=init, interval=50, blit=True
    )

    if save:
        print(f"💾 Guardando animación en {output_name} ...")
        ani.save(output_name, fps=5, dpi=150)
        print("✅ Animación guardada correctamente.")
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Animador de simulación SFM.")
    parser.add_argument("file", help="Ruta al archivo CSV de salida.")
    parser.add_argument("--L", type=float, default=6.0, help="Tamaño del dominio (por defecto 6.0)")
    parser.add_argument("--save", action="store_true", help="Guarda la animación en MP4 en lugar de mostrarla.")
    parser.add_argument("--output", default="animation.mp4", help="Nombre del archivo de salida MP4")

    args = parser.parse_args()
    animate_simulation(args.file, L=args.L, save=args.save, output_name=args.output)
