global using real = System.Double;

namespace nMethods_4;

public enum Test
{
    firstTest,  // две окружности, которые не пересекаются
    secondTest, // две окружности, которые пересекаются в одной точке
    thirdTest,  // две окружности, которые пересекаются в двух точках
    fourthTest, // две окружности, которые пересекаются в двух точках + уравнение прямой
    fifthTest,  // три попарно пересекающиеся прямые
    sixthTest   // прямая, которая пересекает синусоиду 
}

public enum Derivative // метод подсчета производных
{
    Analytic,
    Numerical
}

public class SNE
{
    private delegate void CalcElementsJacobi();
    CalcElementsJacobi calcElementsJacobi;

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
    private Derivative _derivative; // метод подсчета производных
    private real h; // шаг для численного дифференцирования

    // буфер для метода симметризации
    private real[,] temp;
    private real[] tempF;

    public SNE(string path, Test test, Derivative derivative)
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

                x = new Vector(m);
                x.vec = sr.ReadLine().Split(" ").
                        Select(value => real.Parse(value)).ToArray();
            }

            A = new real[m, n];
            temp = new real[n, n];
            tempF = new real[n];
            f = new(m);
            delta = new(m);

            _test = test;
            _derivative = derivative;
            h = 1E-12;

            if (_derivative == Derivative.Analytic)
                calcElementsJacobi = CalcAnalyticElementsJacobi;
            else if (_derivative == Derivative.Numerical)
                calcElementsJacobi = CalcNumericalElementsJacobi;
            else
                throw new ArgumentException(message: "Invalid derivative calculation method",
                                            paramName: nameof(_derivative));
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

            case Test.secondTest:
                A[0, 0] = 2 * (x.vec[0] + 3);
                A[0, 1] = 2 * x.vec[1];
                A[1, 0] = 2 * (x.vec[0] - 1);
                A[1, 1] = 2 * x.vec[1];
                break;

            case Test.thirdTest:
                A[0, 0] = 2 * x.vec[0];
                A[0, 1] = 2 * (x.vec[1] + 2);
                A[1, 0] = 2 * x.vec[0];
                A[1, 1] = 2 * (x.vec[1] - 1);
                break;

            case Test.fourthTest:
                A = ResizeArray(A, m, n);
                A[0, 0] = 2 * x.vec[0];
                A[0, 1] = 2 * (x.vec[1] + 2);
                A[1, 0] = 2 * x.vec[0];
                A[1, 1] = 2 * (x.vec[1] - 1);
                A[2, 0] = 500.0 / 441;
                A[2, 1] = -1;
                break;

            case Test.fifthTest:
                A = ResizeArray(A, m, n);
                A[0, 0] = 1;
                A[0, 1] = -1;
                A[1, 0] = 0.1;
                A[1, 1] = -1;
                A[2, 0] = -1;
                A[2, 1] = -1;
                break;

            case Test.sixthTest:
                A[0, 0] = -8 * Math.Sin(2 * x.vec[0] - 1);
                A[0, 1] = -1;
                A[1, 0] = 1;
                A[1, 1] = -1;
                break;

            default:
                throw new ArgumentException(message: "Invalid enum value",
                                            paramName: nameof(_test));
        }
    }

    private void CalcNumericalElementsJacobi()
    {
        if (m > n)
            A = ResizeArray(A, m, n);

        for (uint i = 0; i < m; i++)
            for (uint j = 0; j < n; j++)
                A[i, j] = Differentiation(i, j);
    }

    private void CalcElementsF()
    {
        for (uint i = 0; i < m; i++)
            f.vec[i] = ValueAtPoint(i, 0, 0);
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

        CalcElementsF();
        primaryNorm = CalcNorm(f.vec);
        currentNorm = previousNorm = CalcNorm(f.vec);

        for (index = 0; index < maxIter &&
        (currentNorm / primaryNorm) >= sndEps; index++)
        {
            previousNorm = CalcNorm(f.vec);

            f = -f;
            beta = 1;

            calcElementsJacobi();

            if (m > n)
                SymmetrizationMethod();

            MethodGauss();

            while (beta >= fstEps)
            {
                x = x + beta * delta;

                CalcElementsF();
                currentNorm = CalcNorm(f.vec);

                if (currentNorm > previousNorm)
                    beta /= 2;
                else
                    break;
            }
        }
    }
    private void SymmetrizationMethod()
    {
        for (uint i = 0; i < n; i++)
            for (uint j = 0; j < m; j++)
                tempF[i] += A[j, i] * f.vec[j];

        for (uint i = 0; i < n; i++)
            for (uint j = 0; j < n; j++)
                for (uint k = 0; k < m; k++)
                    temp[i, j] += A[k, i] * A[k, j];

        A = ResizeArray(A, n, n);
        Copy();
        Array.Copy(tempF, f.vec, n);
        Array.Clear(tempF, 0, n);

        for (int i = 0; i <= n; i++)
            Array.Clear(temp, i, n);
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
                    throw new Exception("Zero element of the column");

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
            for (uint i = 0; i < n; i++)
                sw.WriteLine(x.vec[i]);
        }
    }

    private real ValueAtPoint(uint numberFunc, real hx, real hy)
    {
        return _test switch
        {
            Test.firstTest => numberFunc switch
            {
                0 => (x.vec[0] + hx + 3) * (x.vec[0] + hx + 3) +
                     (x.vec[1] + hy) * (x.vec[1] + hy) - 8,

                1 => (x.vec[0] + hx - 2) * (x.vec[0] + hx - 2) +
                     (x.vec[1] + hy) * (x.vec[1] + hy) - 8,

                _ => throw new ArgumentException(message: "Invalid number function",
                                                 paramName: nameof(numberFunc))
            },

            Test.secondTest => numberFunc switch
            {
                0 => (x.vec[0] + hx + 3) * (x.vec[0] + hx + 3) +
                     (x.vec[1] + hy) * (x.vec[1] + hy) - 4,

                1 => (x.vec[0] + hx - 1) * (x.vec[0] + hx - 1) +
                     (x.vec[1] + hy) * (x.vec[1] + hy) - 4,

                _ => throw new ArgumentException(message: "Invalid number function",
                                                 paramName: nameof(numberFunc))
            },

            Test.thirdTest => numberFunc switch
            {
                0 => (x.vec[0] + hx) * (x.vec[0] + hx) +
                     (x.vec[1] + hy + 2) * (x.vec[1] + hy + 2) - 4,

                1 => (x.vec[0] + hx) * (x.vec[0] + hx) +
                     (x.vec[1] + hy - 1) * (x.vec[1] + hy - 1) - 4,

                _ => throw new ArgumentException(message: "Invalid number function",
                                                 paramName: nameof(numberFunc))
            },

            Test.fourthTest => numberFunc switch
            {
                0 => (x.vec[0] + hx) * (x.vec[0] + hx) +
                     (x.vec[1] + hy + 2) * (x.vec[1] + hy + 2) - 4,

                1 => (x.vec[0] + hx) * (x.vec[0] + hx) +
                     (x.vec[1] + hy - 1) * (x.vec[1] + hy - 1) - 4,

                2 => (500.0 / 441 * x.vec[0] + hx) - (x.vec[1] + hy),

                _ => throw new ArgumentException(message: "Invalid number function",
                                                 paramName: nameof(numberFunc))
            },

            Test.fifthTest => numberFunc switch
            {
                0 => (x.vec[0] + hx) + 1 - (x.vec[1] + hy),

                1 => (0.1 * x.vec[0] + hx) - (x.vec[1] + hy),

                2 => (-x.vec[0] + hx) + 2 - (x.vec[1] + hy),

                _ => throw new ArgumentException(message: "Invalid number function",
                                                 paramName: nameof(numberFunc))
            },

            Test.sixthTest => numberFunc switch
            {
                0 => 2 + 4 * Math.Cos(2 * x.vec[0] + hx + 1) - (x.vec[1] + hy),

                1 => (x.vec[0] + hx) - (x.vec[1] + hy),

                _ => throw new ArgumentException(message: "Invalid number function",
                                                 paramName: nameof(numberFunc))
            },

            _ => throw new ArgumentException(message: "Invalid enum value",
                                             paramName: nameof(_test))
        };
    }

    private real Differentiation(uint numberFunc, uint numberVariable)
    {
        return numberFunc switch
        {
            0 => numberVariable switch
            {
                0 => (ValueAtPoint(numberFunc, h, 0) -
                     ValueAtPoint(numberFunc, -h, 0)) / (2 * h),

                1 => (ValueAtPoint(numberFunc, 0, h) -
                     ValueAtPoint(numberFunc, 0, -h)) / (2 * h),

                _ => throw new ArgumentException(message: "Invalid number variable",
                                                 paramName: nameof(numberVariable))
            },

            1 => numberVariable switch
            {
                0 => (ValueAtPoint(numberFunc, h, 0) -
                     ValueAtPoint(numberFunc, -h, 0)) / (2 * h),

                1 => (ValueAtPoint(numberFunc, 0, h) -
                     ValueAtPoint(numberFunc, 0, -h)) / (2 * h),

                _ => throw new ArgumentException(message: "Invalid number variable",
                                                 paramName: nameof(numberVariable))
            },

            2 => numberVariable switch
            {
                0 => (ValueAtPoint(numberFunc, h, 0) -
                     ValueAtPoint(numberFunc, -h, 0)) / (2 * h),

                1 => (ValueAtPoint(numberFunc, 0, h) -
                     ValueAtPoint(numberFunc, 0, -h)) / (2 * h),

                _ => throw new ArgumentException(message: "Invalid number variable",
                                                 paramName: nameof(numberVariable))
            },

            _ => throw new ArgumentException(message: "Invalid number function",
                                             paramName: nameof(numberFunc))
        };
    }

    private real[,] ResizeArray(real[,] original, int rows, int cols)
    {
        real[,] newArray = new real[rows, cols];
        int minRows = Math.Min(rows, original.GetLength(0));
        int minCols = Math.Min(cols, original.GetLength(1));

        for (int i = 0; i < minRows; i++)
            for (int j = 0; j < minCols; j++)
                newArray[i, j] = original[i, j];

        return newArray;
    }

    private void Copy()
    {
        for (uint i = 0; i < n; i++)
            for (uint j = 0; j < n; j++)
                A[i, j] = temp[i, j];
    }
}