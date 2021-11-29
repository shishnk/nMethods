using System;
using System.IO;
using System.Linq;

using real = System.Double;

namespace nMethods_2
{
    class SLAE
    {
        private uint n; // matrix dimension
        private uint m; // number of zero diagonals
        private real[] mainDiag; // main diagonal
        private real[] l1, l2, l3; // 3 lower diagonals
        private real[] u1, u2, u3; // 3 upper diagonals
        private real[] f; // right-hand vector
        private real[] xk; // primary approximation
        private real[] xk1; // vector for k+1 iteration
        private real[] dif; // xk - x* (where x* is precise solution) 
        private real eps; // solution accuracy 
        private uint countIter; // count of iterations
        private uint maxIter; // maximum number of iterations
        private real conditioningNumber; // number of conditions
        private real normF; // norm of right-hand vector
        private real[] r; // residual
        private real[] product;

        public SLAE(string pathMatrix, string pathVector, string pathParametrs, string pathXk)
        {
            try
            {
                using (var sr = new StreamReader(pathMatrix))
                {
                    // read the matrix parametrs(n and m)
                    uint[] parametrs = sr.ReadLine().Split(" ").Select(value => uint.Parse(value)).ToArray();
                    n = parametrs[0];
                    m = parametrs[1];

                    // read elements of main diagonal
                    mainDiag = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();

                    // read elements of lower diagonals
                    l1 = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();
                    l2 = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();
                    l3 = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();

                    // read elements of upper diagonals
                    u1 = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();
                    u2 = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();
                    u3 = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();
                }

                using (var sr = new StreamReader(pathVector))
                {   // read right-hand vector
                    f = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();
                }

                using (var sr = new StreamReader(pathParametrs))
                {   // read solution accuracy and maximum count of iterations
                    eps = real.Parse(sr.ReadLine());
                    maxIter = uint.Parse(sr.ReadLine());
                }

                using (var sr = new StreamReader(pathXk))
                {   // read vector of primary approximation
                    xk = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();
                }

                xk1 = new real[n];
                dif = new real[n];
                r = new real[n];
                product = new real[n];

                normF = CalcNorm(f); // calculate norm of right-hand vector
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error: {ex.Message}");
            }
        }

        private void Multiply() // Calculate multiply matrix on xk
        {
            for (uint i = 0; i < n; i++)
                product[i] = MultLine(i, xk, 0);
        }

        private real MultLine(uint i, real[] vector, byte method) // multiply element of line on element by vector; 
                                                                  //the third parameter is needed to calculate different sums in the Gauss-Seidel method
        {
            real sum = 0;

            if (method == 0 || method == 1)
            {
                if (i > 0)
                {
                    sum += l1[i - 1] * vector[i - 1];

                    if (i > m + 1)
                        sum += l2[i - m - 2] * vector[i - m - 2];

                    if (i > m + 2)
                        sum += l3[i - m - 3] * vector[i - m - 3];
                }
            }

            if (method == 0 || method == 2)
            {
                sum += mainDiag[i] * vector[i];

                if (i < n - 1)
                {
                    sum += u1[i] * vector[i + 1];

                    if (i < n - m - 2)
                        sum += u2[i] * vector[i + m + 2];

                    if (i < n - m - 3)
                        sum += u3[i] * vector[i + m + 3];
                }
            }

            return sum;
        }

        private real CalcNorm(real[] vector)
        {
            real sum = 0;

            for (uint i = 0; i < vector.Length; i++)
            {
                sum += vector[i] * vector[i];
            }

            return Math.Sqrt(sum);
        }

        public void CalcConditioningNumberEstimation()
        {
            Multiply();

            for (uint i = 0; i < n; i++)
                dif[i] = xk[i] - (i + 1);

            for (uint i = 0; i < n; i++)
                product[i] -= f[i];

            conditioningNumber = (CalcNorm(dif) / (5 * Math.Sqrt(26))) / (CalcNorm(product) / normF);
        }

        private void Jacobi(real w) // Method Jacobi
        {
            for (uint i = 0; i < n; i++)
            {
                real sum = MultLine(i, xk, 0);

                r[i] = f[i] - sum;
                xk1[i] = xk[i] + w * (f[i] - sum) / mainDiag[i];
            }

            Array.Copy(xk1, xk, n);
        }

        private void GaussSeidel(real w) // Method Gauss-Seidel
        {
            for (uint i = 0; i < n; i++)
            {
                real fstSum = MultLine(i, xk1, 1);
                real scdSum = MultLine(i, xk, 2);

                r[i] = f[i] - (fstSum + scdSum);

                xk1[i] = xk[i] + w * r[i] / mainDiag[i];
            }

            Array.Copy(xk1, xk, n);
        }

        public void CalcIter(real w, Method method)
        {
            uint i;

            Array.Clear(xk, 0, xk.Length);

            for (i = 0; i < maxIter; i++)
            {
                Array.Clear(xk1, 0, xk1.Length);

                if (method == 0)
                    Jacobi(w);
                else
                    GaussSeidel(w);

                if (CalcNorm(r) / normF < eps)
                    break;
            }

            countIter = i;
        }

        public void WriteToFile(string path, real w) // Write to file our result
        {
            using (var sw = new StreamWriter(path, true))
            {
                sw.WriteLine($"w = {w:F2} | iterations = {countIter} | number of conditions = {conditioningNumber:F2}");

                for (uint i = 0; i < xk.Length; i++)
                    sw.WriteLine(xk[i]);

                sw.WriteLine("--------------------------------------------------------------");
            }
        }

    }
}