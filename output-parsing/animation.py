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


def animate_simulation(
    file_path,
    L=6.0,
    scale=350,
    save=False,
    maxT=10.0,
    speed=1.0,
    output_dir="animations",
):
    steps = read_simulation_output(file_path)
    n_frames = len(steps)

    # Compute playback duration (T / speed) and interval in ms
    total_playback_time = maxT / speed
    interval = (total_playback_time / n_frames) * 1000  # convert to milliseconds

    print(f"[INFO] Loaded {n_frames} frames.")
    print(f"[INFO] Playback duration: {total_playback_time:.2f} s")
    print(f"[INFO] Frame interval: {interval:.2f} ms")

    fig, ax = plt.subplots()
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_aspect("equal")
    ax.set_title("Pedestrian Simulation (SFM)")

    time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)
    circles = []
    labels = []

    # Initialize circles (one per particle)
    for pid, x, y, vx, vy, r in steps[0][1]:
        color = "red" if pid == 0 else "royalblue"
        c = Circle((x, y), r, color=color, alpha=0.7)
        ax.add_patch(c)
        circles.append(c)
        label = ax.text(
            x,
            y,
            str(pid),
            ha="center",
            va="center",
            fontsize=6,
            color="white",
            weight="bold",
        )
        labels.append(label)

    def init():
        time_text.set_text("")
        return circles + labels + [time_text]

    def update(frame):
        time, particles = steps[frame]
        for c, label, p in zip(circles, labels, particles):
            pid, x, y, vx, vy, r = p
            c.center = (x, y)
            c.radius = r
            label.set_position((x, y))
        time_text.set_text(f"t = {time:.2f} s")
        return circles + labels + [time_text]

    ani = animation.FuncAnimation(
        fig, update, frames=n_frames, init_func=init, interval=interval, blit=True
    )

    if save:
        fps = n_frames / total_playback_time if total_playback_time > 0 else 30
        base_name = os.path.basename(file_path)
        output_path = os.path.join(output_dir, os.path.splitext(base_name)[0] + ".mp4")
        print(f"[INFO] Saving animation to {output_path} (fps={fps:.2f})")
        ani.save(output_path, fps=fps, dpi=150)
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Animate pedestrian simulation output."
    )
    parser.add_argument("file_path", help="Path to simulation output file")
    parser.add_argument("--L", type=float, default=6.0, help="Simulation box size")
    parser.add_argument(
        "--scale", type=float, default=350, help="Animation scaling factor"
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save animation to MP4 instead of showing it",
    )
    parser.add_argument(
        "--maxT", type=float, required=True, help="Total simulated time (s, required)"
    )
    parser.add_argument(
        "--speed",
        type=float,
        default=1.0,
        help="Playback speed factor (1.0 = real-time, 2.0 = double speed, etc.)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="animations",
        help="Output directory",
    )
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
