using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

using real = System.Double;

namespace nMethods_3
{
    class SLAU
    {
        private uint n; // размерность матрицы
        private uint maxIter; // максимальное кол-во итераций
        private uint countIter; // счетчик итераций
        private real eps; // точность решения
        private real squareNorm; // квадрат норм (для условия выхода)
        private Vector<uint> ig; // указатели начала строк
        private Vector<uint> jg; // номера столбцов внедиагональных элементов
        private Vector<real> ggl; // элементы нижнего треугольника
        private Vector<real> ggu; // элементы верхнего треугольника
        private Vector<real> di; // диагональные элементы
        private Vector<real> gglnew; // элементы нижнего при разложении LU
        private Vector<real> ggunew; // элементы верхнего при разложении LU
        private Vector<real> dinew; // диагональные элементы при разложении LU
        private Vector<real> pr; // вектор правой части
        private Vector<real> r; // вектор невязки
        private Vector<real> z; // вектор спуска (сопряженное направление)
        private Vector<real> p; // вспомогательный вектор
        private Vector<real> x; // вектор начального приближения

        public SLAU(string pathParametrs, string pathRows, string pathColumns, string pathLower,
                    string pathUpper, string pathDiag, string pathVector)
        {   // чтение всех данных
            try
            {
                using (var sr = new StreamReader(pathParametrs))
                {   
                    n = uint.Parse(sr.ReadLine());
                    maxIter = uint.Parse(sr.ReadLine());
                    eps = real.Parse(sr.ReadLine());
                }

                using (var sr = new StreamReader(pathRows))
                {
                    ig = new Vector<uint>(n);
                    ig.vec = sr.ReadLine().Split(" ").Select(value => uint.Parse(value)).ToArray();
                }

                using (var sr = new StreamReader(pathColumns))
                {
                    jg = new Vector<uint>(n);
                    jg.vec = sr.ReadLine().Split(" ").Select(value => uint.Parse(value)).ToArray();
                }

                using (var sr = new StreamReader(pathLower))
                {
                    ggl = new Vector<real>(n);
                    ggl.vec = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();
                }

                using (var sr = new StreamReader(pathUpper))
                {
                    ggu = new Vector<real>(n);
                    ggu.vec = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();
                }

                using (var sr = new StreamReader(pathDiag))
                {
                    di = new Vector<real>(n);
                    di.vec = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();
                }

                using (var sr = new StreamReader(pathVector))
                {
                    pr = new Vector<real>(n);
                    pr.vec = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();
                }

                x = new Vector<real>(n);
                r = new Vector<real>(n);
                z = new Vector<real>(n);
                p = new Vector<real>(n);
                gglnew = new Vector<real>((uint)ggl.vec.Length);
                ggunew = new Vector<real>((uint)ggu.vec.Length);
                dinew = new Vector<real>(n);

                Copy(); // копирование элементов для разложения LU

                // используем декремент для прохождения с 0
                for (uint i = 0; i < ig.vec.Length; i++)
                    ig.vec[i]--;

                for (uint i = 0; i < jg.vec.Length; i++)
                    jg.vec[i]--;
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error: {ex.Message}");
            }
        }

        private Vector<real> Mult(Vector<real> vector)
        {
            Vector<real> product = new Vector<real>(n);

            for (uint i = 0; i < n; i++)
            {
                product.vec[i] = di.vec[i] * vector.vec[i];

                for (uint j = ig.vec[i]; j < ig.vec[i + 1]; j++)
                {
                    product.vec[i] += ggl.vec[j] * vector.vec[jg.vec[j]];
                    product.vec[jg.vec[j]] += ggu.vec[j] * vector.vec[i];
                }
            }

            return product;
        }

        private Vector<real> MultDi(Vector<real> vector) // умножение матрицы на вектор 
                                                        // при диагональном предобуславливании,
                                                        // сразу применен ход для решения СЛАУ
        {
            Vector<real> product = new Vector<real>(n);

            for (uint i = 0; i < n; i++)
                product.vec[i] = 1 / Math.Sqrt(di.vec[i]) * vector.vec[i];

            return product;
        }

        public void LOS() // Локально-оптимальная схема без предобуславливания
        {
            uint index;
            real alpha, beta;

            Vector<real> tmp = new Vector<real>(n);

            r = pr - Mult(x);

            Array.Copy(r.vec, z.vec, n);
            Array.Copy(Mult(z).vec, p.vec, n);

            squareNorm = r * r;

            for (index = 0; index < maxIter && squareNorm > eps; index++)
            {
                alpha = (p * r) / (p * p);
                x = x + alpha * z;
                squareNorm = (r * r) - (alpha * alpha) * (p * p);
                r = r - alpha * p;

                tmp = Mult(r);

                beta = -(p * tmp) / (p * p);
                z = r + beta * z;
                p = tmp + beta * p;
            }

            countIter = index;
        }

        public void LOSWithLU() // Локально-оптимальная схема с факторизацией LU
        {
            uint index;
            real alpha, beta;

            Vector<real> tmp = new Vector<real>(n);

            LU();
            r = Direct(pr - MultDi(x));
            z = Reverse(r);
            p = Direct(Mult(z));

            squareNorm = r * r;

            for (index = 0; index < maxIter && squareNorm > eps; index++)
            {
                alpha = (p * r) / (p * p);
                squareNorm = (r * r) - (alpha * alpha) * (p * p);
                x = x + alpha * z;
                r = r - alpha * p;

                tmp = Direct(Mult(Reverse(r)));

                beta = -(p * tmp) / (p * p);
                z = Reverse(r) + beta * z;
                p = tmp + beta * p;
            }

            countIter = index;
        }

        public void LOSWithDi() // Локально-оптимальная схема с диагональным предобуславливанием
        {
            uint index;
            real alpha, beta;

            Vector<real> tmp = new Vector<real>(n);

            r = MultDi(pr - Mult(x));

            z = MultDi(r);
            p = MultDi(Mult(z));

            squareNorm = r * r;

            for (index = 0; index < maxIter && squareNorm > eps; index++)
            {
                alpha = (p * r) / (p * p);
                x = x + alpha * z;
                squareNorm = (r * r) - (alpha * alpha) * (p * p);
                r = r - alpha * p;

                tmp = MultDi(Mult(MultDi(r)));

                beta = -(p * tmp) / (p * p);
                z = MultDi(r) + beta * z;
                p = tmp + beta * p;
            }

            countIter = index;
        }

        private void LU()
        {
            real suml = 0;
            real sumu = 0;
            real sumdi = 0;

            for (uint i = 0; i < n; i++)
            {
                uint i0 = ig.vec[i];
                uint i1 = ig.vec[i + 1];

                for (uint k = i0; k < i1; k++)
                {
                    uint j = jg.vec[k];
                    uint j0 = ig.vec[j];
                    uint j1 = ig.vec[j + 1];
                    uint ik = i0;
                    uint kj = j0;

                    while (ik < k && kj < j1)
                    {
                        if (jg.vec[ik] == jg.vec[kj])
                        {
                            suml += gglnew.vec[ik] * ggunew.vec[kj];
                            sumu += ggunew.vec[ik] * gglnew.vec[kj];
                            ik++;
                            kj++;
                        }
                        else
                        {
                            if (jg.vec[ik] > jg.vec[kj])
                                kj++;
                            else
                                ik++;
                        }
                    }

                    gglnew.vec[k] -= suml;
                    ggunew.vec[k] = (ggunew.vec[k] - sumu) / dinew.vec[j];
                    sumdi += gglnew.vec[k] * ggunew.vec[k];
                    suml = 0;
                    sumu = 0;
                }

                dinew.vec[i] -= sumdi;
                sumdi = 0;
            }
        }

        private Vector<real> Direct(Vector<real> vector)
        {
            Vector<real> y = new Vector<real>(n);

            Array.Copy(vector.vec, y.vec, n);

            real sum = 0;

            for (uint i = 0; i < n; i++)
            {
                uint i0 = ig.vec[i];
                uint i1 = ig.vec[i + 1];

                for (uint k = i0; k < i1; k++)
                    sum += gglnew.vec[k] * y.vec[jg.vec[k]];

                y.vec[i] = (y.vec[i] - sum) / dinew.vec[i];
                sum = 0;
            }

            return y;
        }

        private Vector<real> Reverse(Vector<real> vector)
        {
            Vector<real> result = new Vector<real>(n);

            Array.Copy(vector.vec, result.vec, n);

            for (int i = (int)n - 1; i >= 0; i--)
            {
                uint i0 = ig.vec[i];
                uint i1 = ig.vec[i + 1];

                for (uint k = i0; k < i1; k++)
                    result.vec[jg.vec[k]] -= ggunew.vec[k] * result.vec[i];
            }

            return result;
        }

        public void WriteToFile(string path, string str, real time)
        {
            using (var sw = new StreamWriter(path, true))
            {
                sw.WriteLine($"Method: {str}\n\n iterations: {countIter}, residual: {squareNorm.ToString("0.00E+0")}, time(mcs): {time}\n");

                for (uint i = 0; i < n; i++)
                    sw.WriteLine(x.vec[i]);

                sw.WriteLine("---------------------------------------------------");
            }
        }

        public void Clear()
        {
            Array.Clear(r.vec, 0, r.vec.Length);
            Array.Clear(z.vec, 0, z.vec.Length);
            Array.Clear(p.vec, 0, p.vec.Length);
            Array.Clear(x.vec, 0, x.vec.Length);
        }

        private void Resize()
        {
            Array.Resize<uint>(ref ig.vec, (int)n + 1);
            Array.Resize<real>(ref di.vec, (int)n);
            Array.Resize<real>(ref ggl.vec, (int)(n * (n - 1) / 2));
            Array.Resize<real>(ref ggu.vec, (int)(n * (n - 1) / 2));
            Array.Resize<real>(ref pr.vec, (int)n);
            Array.Resize<real>(ref dinew.vec, (int)n);
            Array.Resize<real>(ref gglnew.vec, (int)(n * (n - 1) / 2));
            Array.Resize<real>(ref ggunew.vec, (int)(n * (n - 1) / 2));
            Array.Resize<real>(ref r.vec, (int)n);
            Array.Resize<real>(ref z.vec, (int)n);
            Array.Resize<real>(ref p.vec, (int)n);
            Array.Resize<real>(ref x.vec, (int)n);
        }

        private void Copy()
        {
            Array.Copy(di.vec, dinew.vec, di.vec.Length);
            Array.Copy(ggl.vec, gglnew.vec, (uint)ggl.vec.Length);
            Array.Copy(ggu.vec, ggunew.vec, (uint)ggu.vec.Length);
        }

        public void GenMatrixHilbert(uint dim)
        {
            n = dim;

            Resize();

            List<uint> jgList = new List<uint>(); // при разной размерности матрицы явно не видно кол-во элементов в jg // TODO 

            ig.vec[0] = 0;

            for (uint i = 0; i < dim; i++)
            {
                di.vec[i] = 1.0 / (2 * (i + 1) - 1);
                ig.vec[i + 1] = ig.vec[i] + i;
            }

            for (uint i = 0; i < dim; i++)
            {
                real sum = 0;
                uint i0 = ig.vec[i];
                uint i1 = ig.vec[i + 1];
                uint j = i - (i1 - i0);

                for (uint k = i0; k < i1; k++, j++)
                {
                    sum = 1.0 / (i + 1 + j + 1 - 1);
                    ggl.vec[k] = sum;
                    ggu.vec[k] = sum;
                    jgList.Insert((int)k, j);
                }
            }

            Array.Resize<uint>(ref jg.vec, jgList.Count());
            jg.vec = jgList.ToArray();
            Copy();
            GenPrHilbert();
        }

        private void GenPrHilbert()
        {
            for (uint i = 0; i < n; i++)
            {
                pr.vec[i] = 0;

                for (real j = 0; j < n; j++)
                    pr.vec[i] += (j + 1.0) / (i + j + 1.0);
            }
        }
    }
}