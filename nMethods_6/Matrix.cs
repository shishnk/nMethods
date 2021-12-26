namespace nMethods_6;

public class Matrix
{
    private double[,] A;
    public int Size { get; init; }

    public double this[int i, int j]
    {
        get => A[i, j];
        set => A[i, j] = value;
    }

    public Matrix(int size)
    {
        Size = size;
        A = new double[size, size];
    }

    public void PrintDense(string path)
    {
        using (var sw = new StreamWriter(path))
        {
            for (int i = 0; i < Size; i++)
            {
                for (int j = 0; j < Size; j++)
                    sw.Write($"{A[i, j].ToString("0.000")}\t\t");

                sw.WriteLine();
            }
        }
    }
    public void LU()
    {
        for (int i = 0; i < Size; i++)
            for (int j = 0; j < Size; j++)
            {
                double suml = 0;
                double sumu = 0;

                if (i < j)
                {
                    for (int k = 0; k < i; k++)
                        sumu += A[i, k] * A[k, j];

                    A[i, j] = (A[i, j] - sumu) / A[i, i];
                }
                else
                {
                    for (int k = 0; k < j; k++)
                        suml += A[i, k] * A[k, j];

                    A[i, j] -= suml;
                }
            }
    }
}