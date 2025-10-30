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

    # Inicializa los círculos (uno por partícula)
    for pid, x, y, vx, vy, r in steps[0][1]:
        color = 'red' if pid == 0 else 'royalblue'
        c = Circle((x, y), r, color=color, alpha=0.7)
        ax.add_patch(c)
        circles.append(c)

    def init():
        time_text.set_text('')
        return circles + [time_text]

    def update(frame):
        time, particles = steps[frame]
        for c, p in zip(circles, particles):
            _, x, y, _, _, r = p
            c.center = (x, y)
            c.radius = r   # las hace un poco más gorditas visualmente
        time_text.set_text(f"t = {time:.2f} s")
        return circles + [time_text]

    ani = animation.FuncAnimation(fig, update, frames=len(steps),
                                  init_func=init, interval=50, blit=True)

    if save:
        ani.save("animation.mp4", fps=30, dpi=150)
    else:
        plt.show()


if __name__ == "__main__":
    animate_simulation("./PedestrianDynamic/outputs/N_10_L6.000/output_N10_L6.000_t20_0000.csv", L=6.0)
