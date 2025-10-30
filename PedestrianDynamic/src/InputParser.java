import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class InputParser {
    private final String inputPath;
    private final int particles_amount;
    public InputParser(String inputPath, int particles_amount) {
        this.inputPath = inputPath;
        this.particles_amount = particles_amount;
    }
    ArrayList<Particle> parseInputs() {
        ArrayList<Particle> particles = new ArrayList<>();
        try {
            List<String> lines = Files.readAllLines(Path.of(inputPath));
            for (String line : lines) {
                String[] parts = line.split(",");
                if (parts.length != 8) {
                    throw new IllegalArgumentException("Invalid input format");
                }
                int id = Integer.parseInt(parts[0]);
                double x = Double.parseDouble(parts[1]);
                double y = Double.parseDouble(parts[2]);
                double v = Double.parseDouble(parts[3]);
                double angle = Double.parseDouble(parts[4]);
                double radius = Double.parseDouble(parts[5]);
                double mass = Double.parseDouble(parts[6]);
                double desiredAngle = Double.parseDouble(parts[7]);
                Particle particle = new Particle(id, x, y, v, angle, radius, mass, desiredAngle);
                particles.add(particle);
            }
            if (particles.size() != particles_amount + 1) { // - 1 por el obstáculo
                throw new IllegalArgumentException("Number of particles does not match the expected amount");
            }
        } catch (Exception e) {
            System.out.println("N: "+ particles_amount + ", lines: "+ particles.size());
            e.printStackTrace();
        }
        return  particles;
    }

}
