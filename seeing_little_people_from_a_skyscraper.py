import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle

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


def animate_simulation(file_path, L=6.0, scale=350, save=False):
    steps = read_simulation_output(file_path)
    fig, ax = plt.subplots()
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_aspect('equal')
    ax.set_title("Pedestrian Simulation (SFM)")

    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

    circles = []
    labels = []

    # Inicializa círculos y textos (uno por partícula)
    for pid, x, y, vx, vy, r in steps[0][1]:
        color = 'red' if pid == 0 else 'royalblue'
        c = Circle((x, y), r, color=color, alpha=0.7)
        ax.add_patch(c)
        circles.append(c)
        label = ax.text(x, y, str(pid),
                        ha='center', va='center',
                        fontsize=6, color='white', weight='bold')
        labels.append(label)

    def init():
        time_text.set_text('')
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

    ani = animation.FuncAnimation(fig, update, frames=len(steps),
                                  init_func=init, interval=50, blit=True)

    ani.save("animation.mp4", fps=5, dpi=150)
    # plt.show()


if __name__ == "__main__":
    animate_simulation("./outputs/N_100_L6.000/output_N100_L6.000_t20_0000.csv", L=6.0)
