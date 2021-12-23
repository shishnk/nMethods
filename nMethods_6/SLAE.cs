namespace nMethods_6;

public static class SLAE
{
    public static Vector<double> Compute(Matrix matrix, Vector<double> vector)
    {
        double max;
        double eps = 1E-14;

        Vector<double> result = new(vector.Length);

        for (int k = 0; k < matrix.Size; k++)
        {
            max = Math.Abs(matrix[k, k]);
            int index = k;

            for (int i = k + 1; i < matrix.Size; i++)
            {
                if (Math.Abs(matrix[i, k]) > max)
                {
                    max = Math.Abs(matrix[i, k]);
                    index = i;
                }
            }

            for (int j = 0; j < matrix.Size; j++)
            {
                (matrix[k, j], matrix[index, j]) =
                    (matrix[index, j], matrix[k, j]);
            }

            (vector[k], vector[index]) = (vector[index], vector[k]);

            for (int i = k; i < matrix.Size; i++)
            {
                double temp = matrix[i, k];

                if (Math.Abs(temp) < eps)
                    throw new Exception("Zero element of the column");

                for (int j = 0; j < matrix.Size; j++)
                    matrix[i, j] /= temp;

                vector[i] /= temp;

                if (i != k)
                {
                    for (int j = 0; j < matrix.Size; j++)
                        matrix[i, j] -= matrix[k, j];

                    vector[i] -= vector[k];
                }
            }
        }

        for (int k = matrix.Size - 1; k >= 0; k--)
        {
            result[k] = vector[k];

            for (int i = 0; i < k; i++)
                vector[i] = vector[i] - matrix[i, k] * result[k];
        }

        return result;
    }
}