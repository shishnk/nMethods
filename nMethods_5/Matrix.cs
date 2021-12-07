global using real = System.Double;

namespace nMethods_5;
public class Matrix
{
    private real[][] A;
    public int Size { get; init; }

    public Matrix(string path)
    {
        try
        {
            using (var sr = new StreamReader(path))
            {
                Size = int.Parse(sr.ReadLine());
                A = sr.ReadToEnd().Split("\n").Select(row => row.Split(" ").Select(value => real.Parse(value)).ToArray()).ToArray();
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine(ex.Message);
        }
    }

    public real this[int i, int j]
    {
        get => A[i][j];
        set => A[i][j] = value;
    }

    public void LU()
    {
        for (int i = 0; i < Size; i++)
        {
            for (int j = 0; j < Size; j++)
            {
                real suml = 0;
                real sumu = 0;

                if (i < j)
                {
                    for (int k = 0; k < i; k++)
                        sumu += A[i][k] * A[k][j];

                    A[i][j] = (A[i][j] - sumu) / A[i][i];
                }
                else
                {
                    for (int k = 0; k < j; k++)
                        suml += A[i][k] * A[k][j];

                    A[i][j] -= suml;
                }
            }
        }
    }
}