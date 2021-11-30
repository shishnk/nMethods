namespace nMethods_5;

public static class SLAE
{
    public static Vector Compute(Matrix matrix, Vector f)
    {
        Vector x = new(f.Length);
        Vector.Copy(f, x);

        for (int i = 0; i < f.Length; i++)
        {
            real sum = 0;

            for (int k = 0; k < i; k++)
                sum += matrix[i, k] * x[k];

            x[i] = (f[i] - sum) / matrix[i, i];
        }

        for (int i = x.Length - 1; i >= 0; i--)
        {
            real sum = 0;

            for (int k = i + 1; k < x.Length; k++)
                sum += matrix[i, k] * x[k];

            x[i] -= sum;
        }

        return x;
    }
}