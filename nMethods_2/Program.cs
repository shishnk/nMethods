namespace nMethods_2
{
    enum Method
    {
        Jacobi,
        GaussSeidel
    }

    class Program
    {
        static void Main(string[] args) // Method CalcIter(first parametr = relaxation value, second parametr = Jacobi || GaussSeidel)
        {
            SLAE slae = new SLAE("matrix.txt", "vector.txt", "parametrs.txt", "xk.txt");

            for (uint w = 1; w <= 121; w += 1)
            {
                slae.CalcIter(w / 100.0, Method.Jacobi);
                slae.CalcConditioningNumberEstimation();
                slae.WriteToFile("result(Jacobi).txt", w / 100.0);
            }

            for (uint w = 1; w < 200; w += 1)
            {
                slae.CalcIter(w / 100.0, Method.GaussSeidel);
                slae.CalcConditioningNumberEstimation();
                slae.WriteToFile("result(GaussSeidel).txt", w / 100.0);
            }
        }
    }
}