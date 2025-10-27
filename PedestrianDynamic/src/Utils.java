import java.util.concurrent.ThreadLocalRandom;

public final class Utils {
    private Utils() {}

    /** [0,1) */
    public static double randDouble01() {
        return ThreadLocalRandom.current().nextDouble();
    }

    /** [min, max) */
    public static double randInRange(double min, double max) {
        return ThreadLocalRandom.current().nextDouble(min, max);
    }

    /** [min, max) (igual que la de C que “excluye derecha”) */
    public static double randInRangeExcludeRight(double min, double max) {
        return randInRange(min, max);
    }

    /** entero uniforme en [min, max] (inclusivo) */
    public static int intRandInRange(int min, int max) {
        return ThreadLocalRandom.current().nextInt(min, max + 1);
    }

    /** módulo “positivo” como en fmod_positive de C */
    public static double fmodPositive(double x, double y) {
        double m = x % y;
        if (m < 0) m += Math.abs(y);
        return m;
    }

}
