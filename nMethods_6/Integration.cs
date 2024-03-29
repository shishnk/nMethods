using static nMethods_6.Spline;

namespace nMethods_6;

public class Integration
{
    Vector<double> points;
    Vector<double> weights;

    public Integration()
    {
        points = new(3);
        weights = new(3);

        points[0] = 0.0;
        points[1] = Math.Sqrt(3.0 / 5);
        points[2] = -Math.Sqrt(3.0 / 5);

        weights[0] = 8.0 / 9;
        weights[1] = 5.0 / 9;
        weights[2] = 5.0 / 9;
    }

    public double GaussOrder5(Basis fstPsi, Basis sndPsi, double x1, double x2)
    {
        double qi, pi;
        double h = 0;
        double result = 0;

        for (int i = 0; i < 3; i++)
        {
            qi = weights[i];
            h = Math.Abs(x2 - x1);
            pi = (x1 + x2 + points[i] * h) / 2.0;

            result += qi * fstPsi(pi, h) * sndPsi(pi, h);
        }

        return result * h / 2.0;
    }
}