namespace nMethods_5;

public class Eigenvalue
{
    private Matrix matrix;
    private Vector vector;
    private real min;
    private real max;
    private int maxIters;
    private int iters;

    public Eigenvalue(Matrix matrix)
    {
        this.matrix = matrix;
        vector = new(matrix.Size);
        vector.Randomize();
        maxIters = 500;

        if (matrix.Size != vector.Length)
            throw new Exception("Matrix and vector of different dimensions");
    }

    public void Compute()
    {
        FindingMaximum();
        FindingMinumum();
    }

    private void FindingMinumum()
    {
        int index;
        real eps = 1E-12;
        real minPrev = 0;

        vector.Randomize();
        matrix.LU();
        Vector nextVector = new(vector.Length);

        for (index = 0; index < maxIters; index++)
        {
            nextVector = SLAE.Compute(matrix, vector);

            if (index % 10 == 0 && index != 0)
                nextVector.Norming();

            minPrev = min;
            min = (nextVector * vector) / (vector * vector);

            Array.Copy(nextVector.vec, vector.vec, nextVector.Length);

            if (Math.Abs(min - minPrev) < eps)
                break;
        }

        min = 1 / min;
        iters = index;
    }

    private void FindingMaximum()
    {
        int index;
        real eps = 1E-12;
        real maxPrev = 0;

        Vector temp = new(vector.Length);

        for (index = 0; ; index++)
        {
            temp = matrix * vector;

            if (index % 10 == 0 && index != 0)
                vector.Norming();

            maxPrev = max;
            max = (temp * vector) / (vector * vector);

            Array.Copy(temp.vec, vector.vec, temp.Length);

            if (Math.Abs(max - maxPrev) < eps)
                break;
        }

        maxIters = index;
    }

    public void WriteToFile(string path)
    {
        using (var sw = new StreamWriter(path))
        {
            sw.WriteLine($"Maximum eigenvalue = {max:F2} | iterations = {maxIters}");
            sw.WriteLine($"Minimum eigenvalue = {min:F2} | iterations = {iters}");
        }
    }
}