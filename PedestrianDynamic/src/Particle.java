import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import static java.lang.Math.*;

public class Particle {
    private final int id;
    private double x;
    private double y;
    private double speed;
    private double angle;
    private final double radius;
    private final double mass;

    private final List<Particle> neighbors;

    public Particle(int id, double x, double y, double speed, double angle, double radius, double mass) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.speed = speed;
        this.angle = angle;
        this.radius = radius;
        this.neighbors = new ArrayList<>();
        this.mass = mass;
    }

    public static Particle cloneOf(Particle p) {
        return new Particle(p.id, p.x, p.y, p.speed, p.angle, p.radius, p.mass);
    }

    public int id() { return id; }
    public double x() { return x; }
    public double y() { return y; }
    public double speed() { return speed; }
    public double angle() { return angle; }
    public double radius() { return radius; }
    public List<Particle> neighbors() { return neighbors; }
    public int neighborCount() { return neighbors.size(); }
    public double getMass() { return mass; }

    public boolean addNeighbor(Particle neighbor) {
        if (neighbor == null || neighbor.id == this.id) return false;
        neighbors.add(neighbor);
        return true;
    }

    private double vx() { return cos(angle) * speed; }
    private double vy() { return sin(angle) * speed; }

    public double distanceEdgeToEdge(Particle other) {
        double dx = abs(x - other.x);
        double dy = abs(y - other.y);
        return sqrt(dx * dx + dy * dy) - radius - other.radius;
    }

    public void update(double boardSize) {
        updatePos(boardSize);
    }


    public void calculateForces() {
    }

    private void updatePos(double boardSize) {
        x = Utils.fmodPositive(x + vx(), boardSize);
        y = Utils.fmodPositive(y + vy(), boardSize);
    }


    /** CSV: id,x,y,vx,vy,angle */
    public String toCsvRow() {
        return String.format(Locale.US,"%d,%f,%f,%f,%f,%f\n", id, x, y, vx(), vy(), radius);
    }

}
