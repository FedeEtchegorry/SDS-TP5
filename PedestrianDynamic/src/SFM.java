import java.util.List;

import static java.lang.Math.*;

public final class SFM {

    private static final double KN = 1.2e5;   // [N/m] fuerza normal (contacto)
    private static final double KT = 2.4e5;   // [kg/m/s] fricción tangencial
    private static final double TAU = 0.5;    // [s] tiempo de relajación
    private static final double VD = 1.7;     // [m/s] velocidad deseada

    public static double[][] force(List<Particle> particles) {

        double[] fx = new double[particles.size()];
        double[] fy = new double[particles.size()];

        for (Particle p1 : particles) {
            for (Particle p2 : p1.neighbors()) {
                if (p1.equals(p2)) continue;
                fx[p1.getId()] += contactForce(particles.get(p1.getId()), particles.get(p2.getId()))[0];
                fy[p1.getId()] += contactForce(particles.get(p1.getId()), particles.get(p2.getId()))[1];
            }
            fx[p1.getId()] += drivingForce(particles.get(p1.getId()))[0];
        }

        return new double[][]{fx, fy};
    }

    private static double[] contactForce(Particle i, Particle j) {
        double dx = j.x() - i.x();
        double dy = j.y() - i.y();
        double dij = hypot(dx, dy);
        if (dij == 0) return new double[]{0, 0};

        double nx = dx / dij;
        double ny = dy / dij;
        double tx = -ny;
        double ty = nx;

        double overlap = i.radius() + j.radius() - dij;
        if (overlap <= 0) return new double[]{0, 0}; // no hay contacto

        double g = overlap; // g(x) = max(0, x)

        // Velocidades
        double vix = cos(i.getAngle()) * i.getSpeed();
        double viy = sin(i.getAngle()) * i.getSpeed();
        double vjx = cos(j.getAngle()) * j.getSpeed();
        double vjy = sin(j.getAngle()) * j.getSpeed();

        double dvx = vjx - vix;
        double dvy = vjy - viy;
        double deltaVt = dvx * tx + dvy * ty;

        // Fuerzas normal y tangencial
        double fn = KN * g;
        double ft = KT * g * deltaVt;

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

    private SFM() {
        throw new IllegalStateException(this.getClass().getName() + " cannot be instantiated.");
    }
}
