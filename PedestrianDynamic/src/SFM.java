import static java.lang.Math.*;

public class SFM {

    private static final double kn = 1.2e5;   // [N/m] fuerza normal (contacto)
    private static final double kt = 2.4e5;   // [kg/m/s] fricción tangencial
    private static final double tau = 0.5;    // [s] tiempo de relajación
    private static final double vd = 1.7;     // [m/s] velocidad deseada

    public static double[] contactForce(Particle i, Particle j) {
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
        double fn = kn * g;
        double ft = kt * g * deltaVt;

        double fx = fn * nx + ft * tx;
        double fy = fn * ny + ft * ty;

        return new double[]{fx, fy};
    }


    public static double[] drivingForce(Particle i) {
        double desiredVx = vd * cos(i.getDesiredAngle());
        double desiredVy = vd * sin(i.getDesiredAngle());

        double vix = cos(i.getAngle()) * i.getSpeed();
        double viy = sin(i.getAngle()) * i.getSpeed();

        double fx = i.getMass() * ((desiredVx - vix) / tau);
        double fy = i.getMass() * ((desiredVy - viy) / tau);

        return new double[]{fx, fy};
    }



}
