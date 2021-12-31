namespace nMethods_6;

public class Spline
{
    public delegate double Basis(double x, double h);
    Basis[] basis, dBasis, ddBasis;
    FiniteElement[] elements;
    Point2D[] points;
    Matrix matrix;
    Vector<double> vector;
    List<Point2D> result;
    Integration integration;
    double alpha;
    double beta;

    public Spline(string pathElements, string pathPoints, string pathParameters)
    {
        try
        {
            using (var sr = new StreamReader(pathElements))
            {
                elements = sr.ReadToEnd().Split("\n").Select(stringElements => FiniteElement.Parse(stringElements)).ToArray();
            }

            using (var sr = new StreamReader(pathPoints))
            {
                points = sr.ReadToEnd().Split("\n").Select(stringPoints => Point2D.Parse(stringPoints)).ToArray();
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

        double x;

        for (int ielem = 0; ielem < elements.Length; ielem++)
            for (int ipoint = 0; ipoint < points.Length; ipoint++)
                if (elements[ielem].Contain(points[ipoint]) && check[ipoint] == 1)
                {
                    check[ipoint] = -1;
                    x = (points[ipoint].X - elements[ielem].LeftBorder) / elements[ielem].Lenght; // $\xi(x) = \dfrac{x - x_i}{h_i}$

                    for (int i = 0; i < basis.Length; i++)
                    {
                        vector[2 * ielem + i] += points[ipoint].Y * basis[i](x, elements[ielem].Lenght);

                        for (int j = 0; j < basis.Length; j++)
                            matrix[2 * ielem + i, 2 * ielem + j] +=
                            basis[i](x, elements[ielem].Lenght) * basis[j](x, elements[ielem].Lenght) +
                            alpha * integration.GaussOrder5(dBasis[i], dBasis[j],
                            elements[ielem].LeftBorder, elements[ielem].RightBorder) +
                            beta * integration.GaussOrder5(ddBasis[i], ddBasis[j],
                            elements[ielem].LeftBorder, elements[ielem].RightBorder);
                    }
                }
    }

    private void ValueAtPoint()
    {
        double x;
        double sum = 0;
        Point2D changed;

        for (int ielem = 0; ielem < elements.Length; ielem++)
        {
            changed = new(elements[ielem].LeftBorder, 0);

            do
            {
                x = (changed.X - elements[ielem].LeftBorder) / elements[ielem].Lenght;

                for (int i = 0; i < basis.Length; i++)
                    sum += vector[2 * ielem + i] * basis[i](x, elements[ielem].Lenght);

                result.Add(new(changed.X, sum));
                changed += (0.2, 0);
                sum = 0;

            } while (elements[ielem].Contain(changed));
        }
    }

    public void WriteToFile(string path)
    {
        using (var sw = new StreamWriter(path))
        {
            for (int i = 0; i < result.Count; i++)
                sw.WriteLine(result[i].ToString());
        }
    }
}