namespace nMethods_5;
public class Vector
{
    private real[] vec;
    public int Length { get; init; }

    public Vector(int dim)
    {
        vec = new real[dim];
        Length = dim;
    }

    public real this[int index]
    {
        get => vec[index];
        set => vec[index] = value;
    }

    public void Randomize()
    {
        for (int i = 0; i < Length; i++)
            vec[i] = new Random().Next(1, 10);
    }

    public static void Copy(Vector source, Vector destination)
    {
        for (int i = 0; i < source.Length; i++)
            destination[i] = source[i];
    }

    public void Norming()
    {
        real norm = CalcNorm();

        for (int i = 0; i < Length; i++)
            vec[i] /= norm;
    }

    private real CalcNorm()
    {
        real result = 0;

        for (int i = 0; i < Length; i++)
            result += vec[i] * vec[i];

        return Math.Sqrt(result);
    }

    public static Vector operator *(Matrix matrix, Vector vector)
    {
        Vector result = new(vector.vec.Length);

        for (int i = 0; i < vector.Length; i++)
            for (int j = 0; j < vector.Length; j++)
                result.vec[i] += matrix[i, j] * vector.vec[j];

        return result;
    }

    public static real operator *(Vector fstVector, Vector sndVector)
    {
        real result = 0;

        for (int i = 0; i < fstVector.Length; i++)
            result += fstVector.vec[i] * sndVector.vec[i];

        return result;
    }
}