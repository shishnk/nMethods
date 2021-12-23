namespace nMethods_6;
public class Spline
{
    public delegate double Basis(double x, double h);

    (int, int)[] elements;
    (double, double)[] points;
    Matrix matrix;
    Vector<double> vector;
    Gauss Gauss;

    public Spline(string pathElements, string pathPoints)
    {
        try
        {
            using (var sr = new StreamReader(pathElements))
            {
                elements = sr.ReadToEnd().Split("\n").Select(row => row.Split(" "))
                           .Select(value => (int.Parse(value[0]), int.Parse(value[1])))
                           .ToArray();
            }

            using (var sr = new StreamReader(pathPoints))
            {
                points = sr.ReadToEnd().Split("\n").Select(row => row.Split(" "))
                           .Select(value => (double.Parse(value[0]), double.Parse(value[1])))
                           .ToArray();
            }

            matrix = new(elements.Length * 2 + 1);
            vector = new(matrix.Size);
            Gauss = new();

            Basis[] basis = new Basis[4]{BasisHermite.Psi1, BasisHermite.Psi2,
                                         BasisHermite.Psi3, BasisHermite.Psi4};
            Basis[] dBasis = new Basis[4]{BasisHermite.dPsi1, BasisHermite.dPsi2,
                                          BasisHermite.dPsi3, BasisHermite.dPsi4};
            Basis[] ddBasis = new Basis[4]{BasisHermite.ddPsi1, BasisHermite.ddPsi2,
                                           BasisHermite.ddPsi3, BasisHermite.ddPsi4};
        }

        catch (Exception ex)
        {
            Console.WriteLine(ex.Message);
        }
    }

    private void Compute()
    {
        AssemblyMatrix();
        AssemblyVector();

    }

    private void AssemblyMatrix()
    {
        double x, h;
    }

    private void AssemblyVector()
    {

    }
}