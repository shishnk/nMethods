global using real = System.Double;

namespace nMethods_4
{
    public class SNE
    {
        private real fstEps; // точность решения СНУ (для беты)
        private real sndEps; // точность решения СНУ (для частного норм)
        private real maxIter; // максимальное кол-во итераций
        private real beta; // параметр итерационного процесса
        private real primaryNorm; // начальная норма вектора правой части
        private Vector x; // начальное приближение
        private Vector f; // вектор правой части
        private Vector delta; // x* - xk (x* - искомое решение)
        private real[,] A; // матрица Якоби

        public SNE(string path)
        {
            try
            { // считывание данных
                using (var sr = new StreamReader(path))
                {
                    fstEps = real.Parse(sr.ReadLine());
                    sndEps = real.Parse(sr.ReadLine());
                    maxIter = real.Parse(sr.ReadLine());

                    x = new Vector(2);
                    x.vec = sr.ReadLine().Split(" ").
                        Select(value => real.Parse(value)).ToArray();
                }

                A = new real[2, 2];
                f = new Vector(2);
                delta = new Vector(2);
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.Message);
            }
        }

        private void CalcElementsJacobi()
        {
            A[0, 0] = 1;
            A[0, 1] = 1;
            A[1, 0] = 2 * x.vec[0];
            A[1, 1] = 2 * x.vec[1];
        }

        private void CalcElementsF()
        {
            f.vec[0] = x.vec[0] + x.vec[1] - 3;
            f.vec[1] = x.vec[0] * x.vec[0] + x.vec[1] * x.vec[1] - 9;
        }

        private real CalcNorm(real[] vector)
        {
            real result = 0;

            for (uint i = 0; i < vector.Length; i++)
                result += vector[i] * vector[i];

            return Math.Sqrt(result);
        }

        public void MethodNewton()
        {
            uint index;
            real currentNorm;
            real previousNorm;
            beta = 1;

            CalcElementsJacobi();
            CalcElementsF();
            primaryNorm = CalcNorm(f.vec);
            currentNorm = previousNorm = CalcNorm(f.vec);

            for (index = 0; index < maxIter && beta >= fstEps &&
                (CalcNorm(f.vec) / primaryNorm) >= sndEps; index++)
            {
                f = -f;

                MethodGauss();

                if (currentNorm > previousNorm)
                    beta /= 2;

                x = x + beta * delta;

                previousNorm = CalcNorm(f.vec);
                CalcElementsF();
                CalcElementsJacobi();
                currentNorm = CalcNorm(f.vec);
            }
        }

        private void MethodGauss()
        {
            real max;
            real eps = 1E-14;

            for (uint k = 0; k < A.GetLength(0); k++)
            {
                max = Math.Abs(A[k, k]);
                uint index = k;

                for (uint i = k + 1; i < A.GetLength(0); i++)
                {
                    if (Math.Abs(A[i, k]) > max)
                    {
                        max = Math.Abs(A[i, k]);
                        index = i;
                    }
                }

                for (uint j = 0; j < A.GetLength(0); j++)
                {
                    (A[k, j], A[index, j]) =
                        (A[index, j], A[k, j]);
                }

                (f.vec[k], f.vec[index]) = (f.vec[index], f.vec[k]);

                for (uint i = k; i < A.GetLength(0); i++)
                {
                    real temp = A[i, k];

                    if (Math.Abs(temp) < eps)
                        throw new Exception("Нулевой элемент столбца");

                    for (uint j = 0; j < A.GetLength(0); j++)
                        A[i, j] /= temp;

                    f.vec[i] /= temp;

                    if (i != k)
                    {
                        for (uint j = 0; j < A.GetLength(0); j++)
                        {
                            A[i, j] -= A[k, j];
                        }

                        f.vec[i] -= f.vec[k];
                    }
                }
            }

            for (int k = A.GetLength(0) - 1; k >= 0; k--)
            {
                delta.vec[k] = f.vec[k];

                for (uint i = 0; i < k; i++)
                    f.vec[i] = f.vec[i] - A[i, k] * delta.vec[k];
            }
        }

        public void WriteToFile(string path)
        {
            using (var sw = new StreamWriter(path))
            {
                for (uint i = 0; i < x.vec.Length; i++)
                    sw.WriteLine(x.vec[i]);
            }
        }
    }
}