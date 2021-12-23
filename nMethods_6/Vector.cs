namespace nMethods_6;

public class Vector<T> where T : INumber<T>
{
    private T[] vec;
    public int Length { get; init; }

    public T this[int index]
    {
        get => vec[index];
        set => vec[index] = value;
    }

    public Vector(int dim)
    {
        vec = new T[dim];
        Length = dim;
    }

    public static T operator *(Vector<T> fstVec, Vector<T> scdVec)
    {
        T result = T.Zero;

        for (int i = 0; i < fstVec.Length; i++)
            result += fstVec.vec[i] * scdVec.vec[i];

        return result;
    }

    public static Vector<T> operator *(double constant, Vector<T> vector)
    {
        Vector<T> result = new(vector.Length);

        for (int i = 0; i < vector.Length; i++)
            result.vec[i] = vector.vec[i] * T.Create(constant);

        return result;
    }

    public static Vector<T> operator +(Vector<T> fstVec, Vector<T> scdVec)
    {
        Vector<T> result = new(fstVec.Length);

        for (int i = 0; i < fstVec.Length; i++)
            result.vec[i] = fstVec.vec[i] + scdVec.vec[i];

        return result;
    }

    public static Vector<T> operator -(Vector<T> fstVec, Vector<T> scdVec)
    {
        Vector<T> result = new(fstVec.Length);

        for (int i = 0; i < fstVec.Length; i++)
            result.vec[i] = fstVec.vec[i] - scdVec.vec[i];

        return result;
    }

    public static void Copy(Vector<T> source, Vector<T> destination)
    {
        for (int i = 0; i < source.Length; i++)
            destination[i] = source[i];
    }
}