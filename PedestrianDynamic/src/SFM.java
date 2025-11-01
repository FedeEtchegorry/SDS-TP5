import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import static java.lang.Math.*;

public final class SFM {

    private static final double KN = 1.2e5;
    private static final double KT = 2.4e5;
    private static final double TAU = 0.5;
    private static final double VD = 1.7;

    private static final ForkJoinPool POOL =
            new ForkJoinPool(Math.max(1, Runtime.getRuntime().availableProcessors()));

    private SFM() { throw new IllegalStateException("nope"); }

    public static double[][] force(List<Particle> particles, double L) {
        final int n = particles.size();
        final int P = Math.max(1, POOL.getParallelism());
        final int STRIPES = Integer.highestOneBit(P * 4);
        final double[][] fxStripe = new double[STRIPES][n];
        final double[][] fyStripe = new double[STRIPES][n];

        try {
            POOL.submit(() ->
                    IntStream.range(1, n).parallel().forEach(i -> {
                        final int stripe = stripeIndex();
                        final double[] fxLoc = fxStripe[stripe];
                        final double[] fyLoc = fyStripe[stripe];

                        final Particle p1 = particles.get(i);

                        double fxLocal = 0.0;
                        double fyLocal = 0.0;

                        for (Particle p2 : p1.neighbors()) {
                            if (p1.notObstacle() && p2.notObstacle() && p1.compareTo(p2) >= 0) continue;

                            double[] c = contactForce(p1, p2, L);
                            fxLocal += c[0];
                            fyLocal += c[1];

                            if (p2.notObstacle()) {
                                int j = p2.getId();
                                fxLoc[j] -= c[0];
                                fyLoc[j] -= c[1];
                            }
                        }

                        double[] d = drivingForce(p1);
                        fxLocal += d[0];
                        fyLocal += d[1];

                        fxLoc[p1.getId()] += fxLocal;
                        fyLoc[p1.getId()] += fyLocal;
                    })
            ).get();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        double[][] f = new double[n][2];
        for (int s = 0; s < STRIPES; s++) {
            double[] sx = fxStripe[s], sy = fyStripe[s];
            for (int i = 0; i < n; i++) {
                f[i][0] += sx[i];
                f[i][1] += sy[i];
            }
        }
        return f;
    }

    private static int stripeIndex() {
        return (int) (Thread.currentThread().getId() & 0x7fffffff) & 0xFF;
    }

    // --- física igual que la tuya ---
    private static double[] contactForce(Particle i, Particle j, double L) {
        double dx = j.x() - i.x();
        double dy = j.y() - i.y();

        if (abs(dx) > L / 2) dx -= signum(dx) * L;
        if (abs(dy) > L / 2) dy -= signum(dy) * L;

        double dij = sqrt(dx * dx + dy * dy);

        double nx = dx / dij, ny = dy / dij;
        double tx = -ny, ty = nx;

        double overlap = i.radius() + j.radius() - dij;
        if (overlap <= 0) return new double[]{0, 0};

        double vix = cos(i.getAngle()) * i.getSpeed();
        double viy = sin(i.getAngle()) * i.getSpeed();
        double vjx = cos(j.getAngle()) * j.getSpeed();
        double vjy = sin(j.getAngle()) * j.getSpeed();

        double dvx = vjx - vix, dvy = vjy - viy;
        double deltaVt = dvx * tx + dvy * ty;

        double fn = -KN * overlap;
        double ft = KT * overlap * deltaVt;

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
