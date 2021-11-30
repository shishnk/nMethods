namespace nMethods_5;

public static class SLAE
{
    public static Vector Compute(Matrix matrix, Vector f)
    {
        Vector x = new(f.Length);
        Array.Copy(f.vec, x.vec, f.Length);

        for (int i = 0; i < f.Length; i++)
        {
            real sum = 0;

            for (int k = 0; k < i; k++)
                sum += matrix.A[i, k] * x.vec[k];

            x.vec[i] = (f.vec[i] - sum) / matrix.A[i, i];
        }

        for (int i = x.Length - 1; i >= 0; i--)
        {
            real sum = 0;

            for (int k = i + 1; k < x.Length; k++)
                sum += matrix.A[i, k] * x.vec[k];

            x.vec[i] -= sum;
        }

        return x;
    }
}