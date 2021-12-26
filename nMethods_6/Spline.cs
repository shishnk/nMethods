namespace nMethods_6;
public class Spline
{
    public delegate double Basis(double x, double h);
    Basis[] basis, dBasis, ddBasis;

    (double, double)[] elements;
    (double, double)[] points;
    Matrix matrix;
    Vector<double> vector;
    Integration integration;
    List<double> result;
    double alpha;
    double beta;


    public Spline(string pathElements, string pathPoints, string pathParameters)
    {
        try
        {
            using (var sr = new StreamReader(pathElements))
            {
                elements = sr.ReadToEnd().Split("\n").Select(row => row.Split(" "))
                           .Select(value => (double.Parse(value[0]), double.Parse(value[1])))
                           .ToArray();
            }

            using (var sr = new StreamReader(pathPoints))
            {
                points = sr.ReadToEnd().Split("\n").Select(row => row.Split(" "))
                           .Select(value => (double.Parse(value[0]), double.Parse(value[1])))
                           .ToArray();
            }

            using (var sr = new StreamReader(pathParameters))
            {
                alpha = double.Parse(sr.ReadLine());
                beta = double.Parse(sr.ReadLine());
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine(ex.Message);
        }

        matrix = new(elements.Length * 2 + 2);
        vector = new(matrix.Size);
        integration = new();
        result = new();

        basis = new Basis[4]{BasisHermite.Psi1, BasisHermite.Psi2,
                                         BasisHermite.Psi3, BasisHermite.Psi4};
        dBasis = new Basis[4]{BasisHermite.dPsi1, BasisHermite.dPsi2,
                                          BasisHermite.dPsi3, BasisHermite.dPsi4};
        ddBasis = new Basis[4]{BasisHermite.ddPsi1, BasisHermite.ddPsi2,
                                           BasisHermite.ddPsi3, BasisHermite.ddPsi4};
    }

    public void Compute()
    {
        Assembly();
        ChangingFunctionality();
        matrix.PrintDense("matrix.txt");
        matrix.LU();

        try
        {
            vector = SLAE.Compute(matrix, vector);
        }
        catch (Exception ex)
        {
            Console.WriteLine(ex.Message);
        }

        ValueAtPoint();
    }

    private void Assembly()
    {
        Vector<int> check = new(points.Length);
        check.Fill(1);

        double x, h;

        for (int ielem = 0; ielem < elements.Length; ielem++)
        {
            h = Math.Abs(elements[ielem].Item2 - elements[ielem].Item1);

            for (int ipoint = 0; ipoint < points.Length; ipoint++)
                if (points[ipoint].Item1 >= elements[ielem].Item1 &&
                    points[ipoint].Item1 <= elements[ielem].Item2 && check[ipoint] == 1)
                {
                    check[ipoint] = -1;

                    x = (points[ipoint].Item1 - elements[ielem].Item1) / h;

                    for (int i = 0; i < basis.Length; i++)
                    {
                        vector[2 * ielem + i] += points[ipoint].Item2 * basis[i](x, h);

                        for (int j = 0; j < basis.Length; j++)
                            matrix[2 * ielem + i, 2 * ielem + j] += basis[i](x, h) * basis[j](x, h);
                    }
                }
        }
    }

    private void ChangingFunctionality() // TODO fix index
    {
        for (int ielem = 0; ielem < elements.Length; ielem++)
            for (int i = 0; i < basis.Length; i++)
                for (int j = 0; j < basis.Length; j++)
                    matrix[i, j] +=
                    alpha * integration.GaussOrder5(dBasis[i], dBasis[j],
                    elements[ielem].Item1, elements[ielem].Item2) +
                    beta * integration.GaussOrder5(ddBasis[i], ddBasis[j],
                    elements[ielem].Item1, elements[ielem].Item2);
    }


    private void ValueAtPoint()
    {
        // List<double> result = new();
        double result = 0;

        double x, h;

        for (int ielem = 0; ielem < elements.Length - 1; ielem++)
        {
            h = Math.Abs(elements[ielem].Item2 - elements[ielem].Item1);
            x = (-1.1 - elements[ielem].Item1) / h;


            for (int i = 0; i < basis.Length; i++)
                result += vector[i] * basis[i](x, h);
        }

        Console.WriteLine(result);
    }

    public void WriteToFile(string path)
    {
        using (var sw = new StreamWriter(path))
        {
            for (int i = 0; i < result.Count; i++)
                sw.WriteLine(result[i]);
        }
    }
}