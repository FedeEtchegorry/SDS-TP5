import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import static java.lang.Math.*;

public final class SFM {

    private static final double KN = 1.2e5;   // [N/m] fuerza normal (contacto)
    private static final double KT = 2.4e5;   // [kg/m/s] fricción tangencial
    private static final double TAU = 0.5;    // [s] tiempo de relajación
    private static final double VD = 1.7;     // [m/s] velocidad deseada

    private SFM() {
        throw new IllegalStateException(this.getClass().getName() + " cannot be instantiated.");
    }

    public static double[][] force(List<Particle> particles, double L) {
        int n = particles.size();
        double[][] f = new double[n][2];

        int availableThreads = Runtime.getRuntime().availableProcessors();
        ForkJoinPool pool = new ForkJoinPool(availableThreads);

        try {
            pool.submit(() ->
                    IntStream.range(1, n).parallel().forEach(i -> {
                        Particle p1 = particles.get(i);

                        double fx = 0, fy = 0;

                        for (Particle p2 : p1.neighbors()) {
                            if (p1.compareTo(p2) >= 0 && p1.notObstacle() && p2.notObstacle()) continue;

                            double[] contact = contactForce(p1, p2, L);

                            fx += contact[0];
                            fy += contact[1];

                            // suma opuesta en p2 (sin tocar el obstáculo id=0)
                            if (p2.notObstacle()) {
                                synchronized (f[p2.getId()]) {
                                    f[p2.getId()][0] -= contact[0];
                                    f[p2.getId()][1] -= contact[1];
                                }
                            }
                        }

                        double[] driving = drivingForce(p1);
                        fx += driving[0];
                        fy += driving[1];

                        f[p1.getId()][0] += fx;
                        f[p1.getId()][1] += fy;
                    })
            ).get();
        } catch (Exception e) {
            throw new RuntimeException(e);
        } finally {
            pool.shutdown();
        }

        return f;
    }

    private static double[] contactForce(Particle i, Particle j, double L) {
        double dx = j.x() - i.x();
        double dy = j.y() - i.y();

        // --- Ajuste por bordes periódicos ---
        if (Math.abs(dx) > L / 2) dx -= Math.signum(dx) * L;
        if (Math.abs(dy) > L / 2) dy -= Math.signum(dy) * L;

        double dij = Math.sqrt(dx * dx + dy * dy);
//        if (dij == 0) return new double[]{0, 0}; // evita división por cero

        double nx = dx / dij;
        double ny = dy / dij;
        double tx = -ny;
        double ty = nx;

        double overlap = i.radius() + j.radius() - dij;
        if (overlap <= 0) return new double[]{0, 0};

        // Velocidades
        double vix = cos(i.getAngle()) * i.getSpeed();
        double viy = sin(i.getAngle()) * i.getSpeed();
        double vjx = cos(j.getAngle()) * j.getSpeed();
        double vjy = sin(j.getAngle()) * j.getSpeed();

        double dvx = vjx - vix;
        double dvy = vjy - viy;
        double deltaVt = dvx * tx + dvy * ty;

        // Fuerzas normal y tangencial
        double fn = -KN * overlap;
        double ft = KT * overlap * deltaVt;

//        if (i.getId() == 0 || j.getId() == 0) {
//            fn *= 50;  // colisión más rígida
//            ft = 0;    // sin fricción
//        }

        double fx = fn * nx + ft * tx;
        double fy = fn * ny + ft * ty;

        return new double[]{fx, fy};
    }


    private static double[] drivingForce(Particle i) {
        double desiredVx = VD * cos(i.getDesiredAngle());
        double desiredVy = VD * sin(i.getDesiredAngle());

        double vix = cos(i.getAngle()) * i.getSpeed();
        double viy = sin(i.getAngle()) * i.getSpeed();

        double fx = i.getMass() * ((desiredVx - vix) / TAU);
        double fy = i.getMass() * ((desiredVy - viy) / TAU);

        return new double[]{fx, fy};
    }
}
