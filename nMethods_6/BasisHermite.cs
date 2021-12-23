namespace nMethods_6;
static public class BasisHermite
{
    public static double Psi1(double x, double h)
    => 1 - 3 * x * x + 2 * x * x * x;

    public static double Psi2(double x, double h)
    => h * (x - 2 * x * x + x * x * x);

    public static double Psi3(double x, double h)
     => 3 * x * x - 2 * x * x * x;

    public static double Psi4(double x, double h)
    => h * (-x * x + x * x * x);

    public static double dPsi1(double x, double h)
    => -6 * (x - x * x) / h;

    public static double dPsi2(double x, double h)
    => 1 - 4 * x + 3 * x * x;

    public static double dPsi3(double x, double h)
    => 6 * (x - x * x) / h;

    public static double dPsi4(double x, double h)
    => -2 * x + 3 * x * x;

    public static double ddPsi1(double x, double h)
    => -6 * (1 - 2 * x) / (h * h);

    public static double ddPsi2(double x, double h)
    => (-4 + 6 * x) / h;

    public static double ddPsi3(double x, double h)
    => 6 * (1 - 2 * x) / (h * h);

    public static double ddPsi4(double x, double h)
    => (-2 + 6 * x) / h;
}