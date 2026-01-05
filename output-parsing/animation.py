import os
import argparse
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle


def read_simulation_output(path):
    steps = []
    with open(path, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        time = float(lines[i].strip())
        i += 1
        particles = []
        while i < len(lines) and "," in lines[i]:
            parts = lines[i].strip().split(",")
            if len(parts) < 6:
                break
            pid = int(parts[0])
            x, y, vx, vy, r = map(float, parts[1:])
            particles.append((pid, x, y, vx, vy, r))
            i += 1
        steps.append((time, particles))
    return steps


def read_collisions(collisions_path):
    """Devuelve un dict {pid: tiempo_primera_colision}."""
    collision_times = {}
    if not os.path.exists(collisions_path):
        print(f"[WARN] No collision file found at {collisions_path}")
        return collision_times

    with open(collisions_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                t, pid = parts
                pid = int(pid)
                if pid == 0:
                    continue
                t = float(t)
                if pid not in collision_times or t < collision_times[pid]:
                    collision_times[pid] = t
    print(f"[INFO] {len(collision_times)} particles collided with the obstacle.")
    return collision_times


def animate_simulation(
        file_path,
        L=12.0,
        scale=350,
        save=False,
        maxT=10.0,
        speed=1.0,
        output_dir="animations",
):
    steps = read_simulation_output(file_path)
    n_frames = len(steps)

    # Archivo de colisiones asociado
    base_name = os.path.splitext(file_path)[0]
    collisions_path = base_name + "_collisions.csv"
    collision_times = read_collisions(collisions_path)

    total_playback_time = maxT / speed
    interval = (total_playback_time / n_frames) * 1000  # ms

    print(f"[INFO] Loaded {n_frames} frames.")
    print(f"[INFO] Playback duration: {total_playback_time:.2f} s")
    print(f"[INFO] Frame interval: {interval:.2f} ms")

    fig, ax = plt.subplots()
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_aspect("equal")
    ax.set_title("Pedestrian Simulation (SFM)")

    time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

    # Diccionario dinámico: pid → (circle, label)
    particles_drawn = {}

    def init():
        time_text.set_text("")
        return [time_text]

    def update(frame):
        time, particles = steps[frame]

        # Actualizar o crear partículas visibles
        current_pids = set()
        for pid, x, y, vx, vy, r in particles:
            current_pids.add(pid)

            if pid not in particles_drawn:
                # Crear nueva partícula
                color = (
                    "red" if pid == 0 else
                    "limegreen" if pid in collision_times and time >= collision_times[pid]
                    else "royalblue"
                )
                c = Circle((x, y), r, color=color, alpha=0.7)
                ax.add_patch(c)
                label = ax.text(
                    x, y, str(pid),
                    ha="center", va="center",
                    fontsize=6, color="white", weight="bold",
                )
                particles_drawn[pid] = (c, label)

            # Actualizar posición/color del pid existente
            c, label = particles_drawn[pid]
            c.center = (x, y)
            c.radius = r
            if pid == 0:
                c.set_color("red")
            elif pid in collision_times and time >= collision_times[pid]:
                c.set_color("limegreen")
            else:
                c.set_color("royalblue")
            label.set_position((x, y))
            label.set_text(str(pid))

        # Eliminar los que ya no están en el frame actual
        missing = set(particles_drawn.keys()) - current_pids
        for pid in missing:
            c, label = particles_drawn[pid]
            c.remove()
            label.remove()
            del particles_drawn[pid]

        time_text.set_text(f"t = {time:.2f} s")

        visible = [obj for pair in particles_drawn.values() for obj in pair]
        return visible + [time_text]

    ani = animation.FuncAnimation(
        fig, update, frames=n_frames, init_func=init,
        interval=interval, blit=True
    )

    if save:
        fps = n_frames / total_playback_time if total_playback_time > 0 else 30
        output_path = os.path.join(output_dir, os.path.basename(base_name) + ".mp4")
        print(f"[INFO] Saving animation to {output_path} (fps={fps:.2f})")
        ani.save(output_path, fps=fps, dpi=150)
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Animate pedestrian simulation output.")
    parser.add_argument("file_path", help="Path to simulation output file")
    parser.add_argument("--L", type=float, default=6.0)
    parser.add_argument("--scale", type=float, default=350)
    parser.add_argument("--save", action="store_true")
    parser.add_argument("--maxT", type=float, required=True)
    parser.add_argument("--speed", type=float, default=1.0)
    parser.add_argument("--output_dir", type=str, default="animations")
    args = parser.parse_args()

    animate_simulation(
        args.file_path,
        L=args.L,
        scale=args.scale,
        save=args.save,
        maxT=args.maxT,
        speed=args.speed,
        output_dir=args.output_dir,
    )
