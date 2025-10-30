import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class ParticlesGenerator {
    private final int particleCount;
    private final double boardSize;
    private final double speed;
    private final double rMin;
    private final double rMax;
    private final double mass;
    private static final String OUTPUT_PATH = "./inputs";

    public ParticlesGenerator(int particleCount, double boardSize, double speed, double rMin, double rMax, double mass) {
        this.particleCount = particleCount;
        this.boardSize = boardSize;
        this.speed = speed;
        this.rMin = rMin;
        this.rMax = rMax;
        this.mass = mass;
    }

    public void generateInputs(int iteration) throws IOException {
        String dirPath = OUTPUT_PATH + String.format("/N%04d", particleCount);
        Files.createDirectories(Path.of(dirPath));
        String fileName = String.format("input_N%04d_%s.csv", particleCount, String.format("%04d", iteration));
        File file = new File(dirPath + "/" + fileName);
        if (file.exists()) {
            file.delete();
        }
        file.createNewFile();

        List<Particle> particles = new ArrayList<>();

    // Obstáculo:
        Particle obstacle = new Particle(0, boardSize / 2, boardSize / 2, 0.0, 0.0, Grid.OBSTACLE_RADIUS, mass, 0.0);
        String center_particle_string = String.format(Locale.US,
                "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.5f,%.17g%n",
                obstacle.id(), obstacle.x(), obstacle.y(), obstacle.getSpeed(), obstacle.getAngle(), obstacle.getRadius(), obstacle.getMass(), obstacle.getDesiredAngle());
        Files.write(file.toPath(), center_particle_string.getBytes(), StandardOpenOption.APPEND);
        particles.add(obstacle);

        // Partículas:
        for (int i = 1; i <= particleCount; i++) {
            double radius = rMin + Math.random() * (rMax - rMin);
            double angle = Math.random() * 2 * Math.PI;
            double desiredAngle = Math.random() * 2 * Math.PI;
            double x, y;
            Particle p;

            int tries = 0;
            boolean valid;
            do {
                x = Math.random() * (boardSize - 2 * radius) + radius;
                y = Math.random() * (boardSize - 2 * radius) + radius;

                p = new Particle(i, x, y, speed, angle, radius, mass, desiredAngle);
                valid = true;

                for (Particle other : particles) {
                    double dx = p.x() - other.x();
                    double dy = p.y() - other.y();
                    double dist2 = dx * dx + dy * dy;
                    double minDist = p.getRadius() + other.getRadius();
                    if (dist2 < minDist * minDist) {
                        valid = false;
                        break;
                    }
                }
                tries++;
                if (tries > 10000) {
                    throw new RuntimeException("No se pudo ubicar la partícula " + i + " sin superposición después de 10000 intentos");
                }
            } while (!valid);

            particles.add(p);

            String particleStr = String.format(Locale.US,
                    "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.5f,%.17g%n",
                    i, x, y, speed, angle, radius, mass, desiredAngle);
            Files.write(file.toPath(), particleStr.getBytes(), StandardOpenOption.APPEND);
        }


        System.out.println("File " + file.getName() + " created successfully.");
    }

    public static void main(String[] args) throws IOException {
        if (args.length < 6) {
            System.out.println("Usage: java ParticlesGenerator <N> <L> <speed> <rMin> <rMax> <iterations> <mass>");
            return;
        }

        int N = Integer.parseInt(args[0]);
        double L = Double.parseDouble(args[1]);
        double speed = Double.parseDouble(args[2]);
        double rMin = Double.parseDouble(args[3]);
        double rMax = Double.parseDouble(args[4]);
        int iterations = Integer.parseInt(args[5]);
        double mass = Double.parseDouble(args[6]);

        if (N <= 0 || L <= 0 || speed <= 0 || rMin <= 0 || rMax <= 0 || iterations <= 0 || rMin > rMax) {
            System.out.println("Error: invalid parameters.");
            return;
        }

        ParticlesGenerator gen = new ParticlesGenerator(N, L, speed, rMin, rMax, mass);
        for (int i = 0; i < iterations; i++) {
            gen.generateInputs(i);
        }
    }
}
