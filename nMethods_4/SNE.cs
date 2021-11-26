global using real = System.Double;

namespace nMethods_4;

public enum Test
{
    firstTest,  // две окружности, которые не пересекаются
    secondTest, // две окружности, которые пересекаются в одной точке
    thirdTest,  // две окружности, которые пересекаются в двух точках
    fourthTest, // три попарно пересекающиеся прямые
    fifthTest,  // исследование влияния взвешивания уравнений СНУ  
    sixthTest   // прямая, которая пересекает синусоиду 
}

public class SNE
{
    private real fstEps; // точность решения СНУ (для беты)
    private real sndEps; // точность решения СНУ (для частного норм)
    private real maxIter; // максимальное кол-во итераций
    private real beta; // параметр итерационного процесса
    private int n; // кол-во переменных уравнений
    private int m;  // кол-во уравнений
    private real primaryNorm; // начальная норма вектора правой части
    private Vector x; // начальное приближение
    private Vector f; // вектор правой части
    private Vector delta; // x* - xk (x* - искомое решение)
    private real[,] A; // матрица Якоби
    private Test _test; // номер теста
    private real h; // шаг для численного дифференцирования

    public SNE(string path, Test test)
    {
        try
        { // считывание и инициализация данных
            using (var sr = new StreamReader(path))
            {
                m = int.Parse(sr.ReadLine());
                n = int.Parse(sr.ReadLine());
                fstEps = real.Parse(sr.ReadLine());
                sndEps = real.Parse(sr.ReadLine());
                maxIter = real.Parse(sr.ReadLine());

                x = new Vector(2);
                x.vec = sr.ReadLine().Split(" ").
                        Select(value => real.Parse(value)).ToArray();
            }

            A = new real[n, n];
            f = new(m);
            delta = new(m);
            _test = test;
            h = 1E-12;
        }
        catch (Exception ex)
        {
            Console.WriteLine(ex.Message);
        }
    }

    private void CalcAnalyticElementsJacobi()
    {
        switch (_test)
        {
            case Test.firstTest:
                A[0, 0] = 2 * (x.vec[0] + 3);
                A[0, 1] = 2 * x.vec[1];
                A[1, 0] = 2 * (x.vec[0] - 2);
                A[1, 1] = 2 * x.vec[1];
                break;
        }
    }

    private void CalcNumericalElementsJacobi()
    {
        switch (_test)
        {
            case Test.firstTest:
                for (uint i = 0; i < n; i++)
                    for (uint j = 0; j < n; j++)
                        A[i, j] = Differentiation(i, j);
                break;
        }
    }

    private void CalcElementsF()
    {
        switch (_test)
        {
            case Test.firstTest:
                for (uint i = 0; i < 2; i++)
                    f.vec[i] = ValueAtPoint(i, 0, 0);
                break;
        }
    }

    private real CalcNorm(real[] vector)
    {
        real result = 0;

        for (uint i = 0; i < vector.Length; i++)
            result += vector[i] * vector[i];

        return Math.Sqrt(result);
    }

    private real ValueAtPoint(uint numberFunc, real hx, real hy)
    {
        // switch (_test)
        // {
        //     case Test.firstTest:
        return numberFunc switch
        {
            0 => (x.vec[0] + hx + 3) * (x.vec[0] + hx + 3) +
                 (x.vec[1] + hy) * (x.vec[1] + hy) - 8,

            1 => (x.vec[0] + hx - 2) * (x.vec[0] + hx - 2) +
                 (x.vec[1] + hy) * (x.vec[1] + hy) - 8,

            _ => throw new Exception("Invalid number")
        };

        // case Test.secondTest:
        //     return numberFunc switch
        //     {
        //         0 => // TODO
        //     };
        // }
    }

    private real Differentiation(uint numberFunc, uint numberVariable)
    {
        switch (_test)
        {
            case Test.firstTest:

                switch (numberFunc)
                {
                    case 0:
                        return numberVariable switch
                        {
                            0 => (ValueAtPoint(numberFunc, h, 0) -
                                 ValueAtPoint(numberFunc, -h, 0)) / (2 * h),

                            1 => (ValueAtPoint(numberFunc, 0, h) -
                                 ValueAtPoint(numberFunc, 0, -h)) / (2 * h),

                            _ => throw new Exception("Invalid number")
                        };

                    case 1:
                        return numberVariable switch
                        {
                            0 => (ValueAtPoint(numberFunc, h, 0) -
                                 ValueAtPoint(numberFunc, -h, 0)) / (2 * h),

                            1 => (ValueAtPoint(numberFunc, 0, h) -
                                 ValueAtPoint(numberFunc, 0, -h)) / (2 * h),

                            _ => throw new Exception("Invalid number")
                        };
                }

                break;
        }
        return 0; // TODO
    }

    public void MethodNewton()
    {
        uint index;
        real currentNorm;
        real previousNorm;
        beta = 1;

        CalcAnalyticElementsJacobi();
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
            CalcAnalyticElementsJacobi();
            currentNorm = CalcNorm(f.vec);
        }
    }

    private void MethodGauss()
    {
        real max;
        real eps = 1E-14;

        for (uint k = 0; k < n; k++)
        {
            max = Math.Abs(A[k, k]);
            uint index = k;

            for (uint i = k + 1; i < n; i++)
            {
                if (Math.Abs(A[i, k]) > max)
                {
                    max = Math.Abs(A[i, k]);
                    index = i;
                }
            }

            for (uint j = 0; j < n; j++)
            {
                (A[k, j], A[index, j]) =
                    (A[index, j], A[k, j]);
            }

            (f.vec[k], f.vec[index]) = (f.vec[index], f.vec[k]);

            for (uint i = k; i < n; i++)
            {
                real temp = A[i, k];

                if (Math.Abs(temp) < eps)
                    throw new Exception("Нулевой элемент столбца");

                for (uint j = 0; j < n; j++)
                    A[i, j] /= temp;

                f.vec[i] /= temp;

                if (i != k)
                {
                    for (uint j = 0; j < n; j++)
                    {
                        A[i, j] -= A[k, j];
                    }

                    f.vec[i] -= f.vec[k];
                }
            }
        }

        for (int k = n - 1; k >= 0; k--)
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