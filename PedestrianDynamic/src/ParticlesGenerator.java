import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
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

    // Obstáculo:
        String center_particle_string = String.format(Locale.US,
                "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.5f%n",
                0, boardSize / 2, boardSize / 2, 0.0, 0.0, Grid.OBSTACLE_RADIUS, Double.MAX_VALUE);
        Files.write(file.toPath(), center_particle_string.getBytes(), StandardOpenOption.APPEND);

    // Partículas:
        for (int i = 1; i <= particleCount; i++) {
            double radius = rMin + Math.random() * (rMax - rMin);
            double x = Math.random() * (boardSize - 2 * radius) + radius;
            double y = Math.random() * (boardSize - 2 * radius) + radius;
            double angle = Math.random() * 2 * Math.PI;
            double vx = speed * Math.cos(angle);
            double vy = speed * Math.sin(angle);

            String particle = String.format(Locale.US,
                    "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.5f%n",
                    i, x, y, speed, angle, radius, mass);
            Files.write(file.toPath(), particle.getBytes(), StandardOpenOption.APPEND);
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
