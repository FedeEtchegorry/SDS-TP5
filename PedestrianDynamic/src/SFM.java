import java.util.List;

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

        double[][] f = new double[particles.size()][2];

        for (Particle p1 : particles.subList(1, particles.size())) {
            for (Particle p2 : p1.neighbors()) {
                if (p1.equals(p2)) continue;

                f[p1.getId()] = contactForce(p1, p2, L);
            }

            double[] drivingF = drivingForce(p1);
            f[p1.getId()][0] += drivingF[0];
            f[p1.getId()][1] += drivingF[1];
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

        if (i.getId() == 0 || j.getId() == 0) {
            fn *= 50;  // colisión más rígida
            ft = 0;    // sin fricción
        }

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
