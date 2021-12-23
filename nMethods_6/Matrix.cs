namespace nMethods_6;

public class Matrix
{
    private double[][] A;
    public int Size { get; init; }

    public double this[int i, int j]
    {
        get => A[i][j];
        set => A[i][j] = value;
    }

    public Matrix(int size)
    {
        Size = size;
    }
}