import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.*;

public class Grid {
    private final List<Particle> particles;
    private final int N;
    private final Cell[] grid;
    private final double cellSize;
    private final double L;
    private final int M;
    private final double Rc;
    private final boolean periodicBorders;

    public Grid(double L, double Rc, int N, List<Particle> particles, boolean periodicBorders) {
        this.particles = particles;
        this.N = N;
        this.L = L;
        this.Rc = Rc;
        this.periodicBorders = periodicBorders;

        double maxR = 0.0;
        for (Particle p : particles) maxR = Math.max(maxR, p.radius());
        int m = (int) (L / (Rc + 2 * maxR));
        if (m < 1) m = 1;

        this.M = m;
        this.grid = new Cell[M * M];
        for (int i = 0; i < grid.length; i++) grid[i] = new Cell();

        this.cellSize = L / M;
    }

    public void findNeighbors() {
        addParticlesToCells();
        for (int r = 0; r < M; r++) {
            for (int c = 0; c < M; c++) {
                findCellNeighbors(r, c);
            }
        }
    }

    public double getL(){
        return L;
    }

    private void addParticlesToCells() {
        for (Cell c : grid) c.clear();
        for (Particle p : particles) addToCell(p);
    }

    private void addToCell(Particle p) {
        int row = (int) (p.y() / cellSize);
        int col = (int) (p.x() / cellSize);
        Cell cell = getCell(row, col);
        if (cell == null) throw new IllegalStateException("Cell is null at (" + row + "," + col + ")");
        synchronized (cell) {
            cell.add(p);
        }
    }

    private void findCellNeighbors(int row, int col) {
        Cell c = getCell(row, col);
        if (c == null) return;

        // vecinos (4 celdas: arriba, arriba-der, der, abajo-der) + la propia
        Cell[] nbs = new Cell[]{
                c,
                getCell(row - 1, col),
                getCell(row - 1, col + 1),
                getCell(row, col + 1),
                getCell(row + 1, col + 1)
        };

        // 1) Todos-vs-todos dentro de la celda propia
        final List<Particle> ps = nbs[0].particles;
        for (int i = 0; i < ps.size(); i++) {
            Particle a = ps.get(i);
            for (int j = i + 1; j < ps.size(); j++) {
                Particle b = ps.get(j);
                if (calculateDistance(a, b) <= Rc) {
                    a.addNeighbor(b);
                    b.addNeighbor(a);
                }
            }
        }

        // 2) Partículas de la celda propia vs cada celda vecina
        for (int idx = 1; idx < nbs.length; idx++) {
            Cell nc = nbs[idx];
            if (nc == null) continue;
            for (Particle a : ps) {
                for (Particle b : nc.particles) {
                    if (calculateDistance(a, b) <= Rc) {
                        a.addNeighbor(b);
                        b.addNeighbor(a);
                    }
                }
            }
        }
    }

    /**
     * Distancia con bordes periódicos + radios (>=0).
     */
    private double calculateDistance(Particle p1, Particle p2) {
        double dx = abs(p1.x() - p2.x());
        double dy = abs(p1.y() - p2.y());

        if (periodicBorders) {
            dx = min(dx, L - dx);
            dy = min(dy, L - dy);
        }
        return max(0.0, sqrt(dx * dx + dy * dy) - p1.radius() - p2.radius());
    }

    private Cell getCell(int row, int col) {
        if (periodicBorders && M > 2) {
            row = (row % M + M) % M;
            col = (col % M + M) % M;
        } else if (row < 0 || col < 0 || row >= M || col >= M) {
            return null;
        }
        return grid[row * M + col];
    }

    private static final class Cell {
        final List<Particle> particles = new ArrayList<>();

        void clear() {
            particles.clear();
        }

        void add(Particle p) {
            particles.add(p);
        }
    }
}
