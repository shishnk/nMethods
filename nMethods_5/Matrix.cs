global using real = System.Double;

namespace nMethods_5;
public class Matrix
{
    public real[,] A;
    public int Size { get; init; }

    public Matrix(string path)
    {
        try
        {
            using (var sr = new StreamReader(path))
            {
                Size = int.Parse(sr.ReadLine());

                string buffer = sr.ReadToEnd();
                string[] row = buffer.Split('\n');
                string[] column = row[0].Split(' ');

                A = new real[Size, Size];

                real element = 0;

                for (int i = 0; i < Size; i++)
                {
                    column = row[i].Split(' ');

                    for (int j = 0; j < Size; j++)
                    {
                        element = real.Parse(column[j]);
                        A[i, j] = element;
                    }
                }
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine(ex.Message);
        }
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
}