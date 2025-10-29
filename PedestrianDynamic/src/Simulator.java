import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Locale;

public class Simulator {

    private final Grid grid;
    private final ArrayList<Particle> particles;
    private final double deltaT;
    private double t=0.0;
    private final double maxT;
    private final double writeInterval;
    private double nextWriteTime=0;
    private final Path outputPath;


    public Simulator(double L, ArrayList<Particle> particleList, double deltaT, Path outputPath, int maxT) throws IOException {
        this.grid = new Grid(L, particleList.getFirst().radius()*2,particleList.size(), particleList, true);
        this.particles = particleList;
        this.deltaT = deltaT;
        this.maxT = maxT;
        this.writeInterval = 0.1;
        this.outputPath=outputPath;
        executeSimulation();

    }
    public void executeSimulation() throws IOException {
        try (OutputWriter out = OutputWriter.open(outputPath)){


            double[][][] gearR = new double[2][][];
            double[][][] gearP = new double[2][][];

            SFMIntegrator.initGear(particles, deltaT, gearR, gearP);

            while (t<maxT)  {
                if (t>=nextWriteTime){
                    out.writeStep(particles, t);
                    nextWriteTime+=writeInterval;
                }
                grid.findNeighbors();
                SFMIntegrator.updateParticlesGear5(particles, deltaT, gearR, gearP, grid.getL());

                t+=deltaT;
            }
        }
    }


    public static void main(String[] args) throws IOException {
        if (args.length != 7) {
            throw new IllegalArgumentException("Wrong number of arguments");
        }
        int N = Integer.parseInt(args[0]);
        double L = Double.parseDouble(args[1]);
        int iterations = Integer.parseInt(args[2]);
        String inputDir = args[3];
        String outputDir = args[4];
        int maxT = Integer.parseInt(args[5]);
        double deltaT = Double.parseDouble(args[6]);
        if (N <= 0 || L <= 0 || iterations <= 0 || maxT <= 0 || deltaT<=0 || inputDir.isEmpty() || outputDir.isEmpty()) {
            System.out.println("Error: Parameters should be: N L iterations inputDir outputDir maxT deltaT");
            return;
        }
        for (int i = 0; i < iterations; i++) {
            InputParser parser = new InputParser(inputDir + "/N" + String.format("%04d", N) +
                    "/input_N" + String.format("%04d", N) + "_" + String.format("%04d", i) + ".csv", N);
            ArrayList<Particle> particles = parser.parseInputs();
            String L_dir = String.format(Locale.US, "L%.3f", L);
            Path directory = Files.createDirectories(Path.of(outputDir, "N_" + N + "_" + L_dir));
            Path fileName = Path.of(directory + String.format("/output_N%d_%s_t%d_%s.csv", N, L_dir, maxT, String.format("%04d", i)));
            System.out.printf("\nStarting iteration %d/%d...\n", i + 1, iterations);
            Simulator s = new Simulator(L, particles, deltaT,fileName, maxT);
            System.out.println("\nIteration " + (i + 1) + " completed.");
        }
    }
}
