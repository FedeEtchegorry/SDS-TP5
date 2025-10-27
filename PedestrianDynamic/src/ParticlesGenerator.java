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
    private final double radius;
    private final double mass;
    private static final String OUTPUT_PATH = "./inputs";

    public ParticlesGenerator(int particleCount, double boardSize, double speed, double radius, double mass) {
        this.particleCount = particleCount;
        this.boardSize = boardSize;
        this.speed = speed;
        this.radius = radius;
        this.mass = mass;
    }

//    public boolean checkOverlap(Particle p1, Particle p2) {
//        double dx = p1.getX() - p2.getX();
//        double dy = p1.getY() - p2.getY();
//        double distance = Math.sqrt(dx * dx + dy * dy);
//        return distance < (p1.getRadius() + p2.getRadius());
//    }
//
//    public boolean checkCollision(Particle[] particles, Particle new_particle, int created_particles) {
//        for (int i = 0; i < created_particles; i++) {
//            if (checkOverlap(particles[i], new_particle)) return true;
//        }
//        return false;
//    }


    public void generateInputs(int iteration) throws IOException {
        String dirPath = OUTPUT_PATH + String.format("/N%04d", particleCount);
        Files.createDirectories(Path.of(dirPath));
        String fileName = String.format("input_N%04d_%s.csv", particleCount, String.format("%04d", iteration));
        File file = new File(dirPath + "/" + fileName);
        if (file.exists()) {
            file.delete();
        }
        file.createNewFile();
        int i = 0;
        while (i < particleCount) {
            double x = Math.random() * (boardSize - 2 * radius) + radius;
            double y = Math.random() * (boardSize - 2 * radius) + radius;
            double angle = Math.random() * 2 * Math.PI;
//            double vx = speed * Math.cos(angle);
//            double vy = speed * Math.sin(angle);
//            Particle new_particle = new Particle(i,x, y, speed, angle, radius, particleCount-1, mass);
            String particle = String.format(Locale.US, "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.5f%n",i, x, y, speed, angle, radius, mass);
            Files.write(file.toPath(), particle.getBytes(), StandardOpenOption.APPEND);
            i++;

        }
        System.out.println("File " + file.getName() + " created successfully.");
    }


    public static void main(String[] args) throws IOException {
        int N = Integer.parseInt(args[0]);
        double L = Double.parseDouble(args[1]);
        double speed = Double.parseDouble(args[2]);
        double radius = Double.parseDouble(args[3]);
        int iterations = Integer.parseInt(args[4]);
        int mass = Integer.parseInt(args[5]);
        if (N <= 0 || L <= 0 || speed <= 0 || radius <= 0 || iterations <= 0) {
            System.out.println("Error: Parameters should be: N, L, speed, radius, iterations, mass");
            return;
        }
        ParticlesGenerator gen = new ParticlesGenerator(N, L, speed, radius, mass);
        for (int i = 0; i < iterations; i++) {
            gen.generateInputs(i);
        }
    }
}

