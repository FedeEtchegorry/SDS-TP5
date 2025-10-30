import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Random;

import static java.lang.Math.*;

public class Particle {
    private final int id;
    private int logicId;
    private double x;
    private double y;
    private double speed;
    private double angle;
    private boolean collisionWithCenterHasHappened;
    private final double radius;
    private final double mass;
    private final double desiredAngle;

    private final List<Particle> neighbors;

    public Particle(int id, double x, double y, double speed, double angle, double radius, double mass, double desiredAngle) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.speed = speed;
        this.angle = angle;
        this.radius = radius;
        this.mass = mass;
        this.neighbors = new ArrayList<>();
        this.desiredAngle = desiredAngle;
        collisionWithCenterHasHappened = false;
        logicId = id;
    }
    public static Particle cloneOf(Particle p) {
        return new Particle(p.id, p.x, p.y, p.speed, p.angle, p.radius, p.mass, p.desiredAngle);
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

    public boolean getCollisionWithCenterHasHappened() { return collisionWithCenterHasHappened; }
    public void setCollisionWithCenterHasHappened() {
        collisionWithCenterHasHappened = true;
    }
    private void resetCollisionWithCenterHasHappened() {
        collisionWithCenterHasHappened = false;
    }
    public void changeId (int newId) {
        logicId = newId;
        resetCollisionWithCenterHasHappened();
    }

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
        return String.format(Locale.US,"%d,%f,%f,%f,%f,%f\n", logicId, x, y, vx(), vy(), radius);
    }


    public List<Particle> getNeighbors() {
        return neighbors;
    }

    @Override
    public boolean equals(Object obj) {
        return obj instanceof Particle && getId() == ((Particle) obj).getId();
    }

    public double getDesiredAngle() {
        return desiredAngle;
    }

    public double getRadius() {
        return radius;
    }

    public double getAngle() {
        return angle;
    }

    public double getSpeed() {
        return speed;
    }

    public double getY() {
        return y;
    }

    public double getX() {
        return x;
    }

    public int getId() {
        return id;
    }

    public void setPosition(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public void setSpeed(double speed) {
        this.speed = speed;
    }

    public void setAngle(double angle) {
        this.angle = angle;
    }

    public void setVelocity(double vx, double vy) {
        setSpeed(sqrt(vx*vx + vy*vy));
        setAngle(atan2(vy, vx));
    }

    public int getLogicId() {
        return logicId;
    }

    public boolean overlapsWithParticle(Particle other) {
        if (other == null || other.equals(this)) return false;
        double dx = this.x - other.x;
        double dy = this.y - other.y;
        double r  = this.radius + other.radius;
        return dx * dx + dy * dy <= r * r;
    }
}
