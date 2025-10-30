import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class SFMIntegrator {

    private static int maxN;

    /** Inicializa los coeficientes del método Gear-5 para todas las partículas */
    public static void initGear(List<Particle> particles, double[][][] r, double[][][] p, double L) {
        int N = particles.size();
        // r[dim][k][i] → dim={0:x,1:y}, k=0..5
        r[0] = new double[6][N];
        r[1] = new double[6][N];
        p[0] = new double[6][N];
        p[1] = new double[6][N];

        // Posiciones y velocidades iniciales
        for (int i = 1; i < N; i++) {
            Particle pi = particles.get(i);
            r[0][0][i] = pi.getX();
            r[1][0][i] = pi.getY();
            double vx = Math.cos(pi.getAngle()) * pi.getSpeed();
            double vy = Math.sin(pi.getAngle()) * pi.getSpeed();
            r[0][1][i] = vx;
            r[1][1][i] = vy;
        }

        // Aceleraciones iniciales = F/m
        double[][] forces = SFM.force(particles, L);
        for (int i = 1; i < N; i++) {
            Particle pi = particles.get(i);
            double fx = forces[i][0];
            double fy = forces[i][1];
            r[0][2][i] = fx / pi.getMass();
            r[1][2][i] = fy / pi.getMass();
        }

        // Derivadas de orden 3–5 se inicializan en 0
        for (int k = 3; k <= 5; k++) {
            Arrays.fill(r[0][k], 0.0);
            Arrays.fill(r[1][k], 0.0);
        }
        maxN = N;
    }

    /** Un paso del integrador Gear-5 (basado en diap. de teórica y TP previos) */
    public static void updateParticlesGear5(List<Particle> particles, double dt,
                                            double[][][] r, double[][][] p, double boardSize) {
        int N = particles.size();
        double dt1 = dt, dt2 = dt1*dt, dt3 = dt2*dt, dt4 = dt3*dt, dt5 = dt4*dt;

        // Coeficientes del Gear
        final double C0 = GearCoefficients.C0;
        final double C1 = GearCoefficients.C1;
        final double C2 = GearCoefficients.C2;
        final double C3 = GearCoefficients.C3;
        final double C4 = GearCoefficients.C4;
        final double C5 = GearCoefficients.C5;

        int[] kxOld = new int[N], kyOld = new int[N];
        for (int i = 1; i < N; i++) {
            kxOld[i] = (int)Math.floor(r[0][0][i] / boardSize);
            kyOld[i] = (int)Math.floor(r[1][0][i] / boardSize);
        }

        // --- PREDICTOR ---
        for (int i = 1; i < N; i++) {
            for (int d = 0; d < 2; d++) { // 0:x, 1:y
                double r0 = r[d][0][i];
                double r1 = r[d][1][i];
                double r2 = r[d][2][i];
                double r3 = r[d][3][i];
                double r4 = r[d][4][i];
                double r5 = r[d][5][i];

                p[d][0][i] = r0 + r1*dt1 + r2*dt2/2 + r3*dt3/6 + r4*dt4/24 + r5*dt5/120;
                p[d][1][i] = r1 + r2*dt1 + r3*dt2/2 + r4*dt3/6 + r5*dt4/24;
                p[d][2][i] = r2 + r3*dt1 + r4*dt2/2 + r5*dt3/6;
                p[d][3][i] = r3 + r4*dt1 + r5*dt2/2;
                p[d][4][i] = r4 + r5*dt1;
                p[d][5][i] = r5;
            }
        }

        // Actualizar posiciones predichas en las partículas
        for (int i = 1; i < N; i++) {
            Particle pi = particles.get(i);
            pi.setPosition(Utils.fmodPositive(p[0][0][i], boardSize),
                           Utils.fmodPositive(p[1][0][i], boardSize));
        }

        // --- EVALUAR NUEVAS FUERZAS ---
        double[][] F = SFM.force(particles, boardSize);
        double[] R2x = new double[N];
        double[] R2y = new double[N];

        for (int i = 1; i < N; i++) {
            Particle pi = particles.get(i);
            double axPred = p[0][2][i];
            double ayPred = p[1][2][i];
            double axReal = F[i][0] / pi.getMass();
            double ayReal = F[i][1] / pi.getMass();
            R2x[i] = 0.5 * (axReal - axPred) * dt2;
            R2y[i] = 0.5 * (ayReal - ayPred) * dt2;
        }

        // --- CORRECTOR ---
        double invDt  = 1.0 / dt;
        double invDt2 = invDt * invDt;
        double invDt3 = invDt2 * invDt;
        double invDt4 = invDt3 * invDt;
        double invDt5 = invDt4 * invDt;

        for (int i = 1; i < N; i++) {
            for (int d = 0; d < 2; d++) {
                double R2 = (d == 0) ? R2x[i] : R2y[i];

                r[d][0][i] = p[d][0][i] + C0 * R2;
                r[d][1][i] = p[d][1][i] + C1 * R2 * invDt;
                r[d][2][i] = p[d][2][i] + C2 * 2.0 * R2 * invDt2;
                r[d][3][i] = p[d][3][i] + C3 * 6.0 * R2 * invDt3;
                r[d][4][i] = p[d][4][i] + C4 * 24.0 * R2 * invDt4;
                r[d][5][i] = p[d][5][i] + C5 * 120.0 * R2 * invDt5;
            }

            int kxNew = (int)Math.floor(r[0][0][i] / boardSize);
            int kyNew = (int)Math.floor(r[1][0][i] / boardSize);
            boolean exited = (kxNew != kxOld[i]) || (kyNew != kyOld[i]);
            if (exited) {
                particles.get(i).changeId(++maxN);
            }

            // Actualizar partícula
            double vx = r[0][1][i];
            double vy = r[1][1][i];
            Particle pi = particles.get(i);
            pi.setPosition(Utils.fmodPositive(r[0][0][i], boardSize),
                           Utils.fmodPositive(r[1][0][i], boardSize));
            pi.setVelocity(vx, vy);
        }
    }

    private SFMIntegrator() {
        throw new IllegalStateException(this.getClass().getSimpleName() + " cannot be instantiated.");
    }
}
