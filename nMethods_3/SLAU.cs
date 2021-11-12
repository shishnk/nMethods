using System;
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
        private real normPr; // норма вектора правой части
        private uint[] ig; // указатели начала строк
        private uint[] jg; // номера столбцов внедиагональных элементов
        private Vector ggl; // элементы нижнего треугольника
        private Vector ggu; // элементы верхнего треугольника
        private Vector di; // диагональные элементы
        private Vector gglnew; // элементы нижнего при разложении LU
        private Vector ggunew; // элементы верхнего при разложении LU
        private Vector dinew; // диагональные элементы при разложении LU
        private Vector pr; // вектор правой части
        private Vector r; // вектор невязки
        private Vector z; // вектор спуска (сопряженное направление)
        private Vector p; // вспомогательный вектор
        private Vector x; // вектор начального приближения

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
                    ig = sr.ReadToEnd().Split("\n").Select(value => uint.Parse(value)).ToArray();
                }

                using (var sr = new StreamReader(pathColumns))
                {
                    jg = sr.ReadToEnd().Split("\n").Select(value => uint.Parse(value)).ToArray();
                }

                using (var sr = new StreamReader(pathLower))
                {
                    ggl = new Vector(n);
                    ggl.vec = sr.ReadToEnd().Split("\n").Select(value => real.Parse(value)).ToArray();
                }

                using (var sr = new StreamReader(pathUpper))
                {
                    ggu = new Vector(n);
                    ggu.vec = sr.ReadToEnd().Split("\n").Select(value => real.Parse(value)).ToArray();
                }

                using (var sr = new StreamReader(pathDiag))
                {
                    di = new Vector(n);
                    di.vec = sr.ReadToEnd().Split("\n").Select(value => real.Parse(value)).ToArray();
                }

                using (var sr = new StreamReader(pathVector))
                {
                    pr = new Vector(n);
                    pr.vec = sr.ReadToEnd().Split("\n").Select(value => real.Parse(value)).ToArray();
                }

                x = new Vector(n);
                r = new Vector(n);
                z = new Vector(n);
                p = new Vector(n);
                gglnew = new Vector((uint)ggl.vec.Length);
                ggunew = new Vector((uint)ggu.vec.Length);
                dinew = new Vector(n);

                Copy(); // копирование элементов для разложения LU, LL^T
                normPr = CalcNorm(pr);

                // используем декремент для прохождения с 0
                for (uint i = 0; i < ig.Length; i++)
                    ig[i]--;

                for (uint i = 0; i < jg.Length; i++)
                    jg[i]--;
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error: {ex.Message}");
            }
        }

        private Vector Mult(Vector vector)
        {
            Vector product = new Vector(n);

            for (uint i = 0; i < n; i++)
            {
                product.vec[i] = di.vec[i] * vector.vec[i];

                for (uint j = ig[i]; j < ig[i + 1]; j++)
                {
                    product.vec[i] += ggl.vec[j] * vector.vec[jg[j]];
                    product.vec[jg[j]] += ggu.vec[j] * vector.vec[i];
                }
            }

            return product;
        }

        private Vector MultT(Vector vector)
        {
            Vector product = new Vector(n);

            for (uint i = 0; i < n; i++)
            {
                product.vec[i] = di.vec[i] * vector.vec[i];

                for (uint j = ig[i]; j < ig[i + 1]; j++)
                {
                    product.vec[i] += ggu.vec[j] * vector.vec[jg[j]];
                    product.vec[jg[j]] += ggl.vec[j] * vector.vec[i];
                }
            }

            return product;
        }

        private Vector MultDi(Vector vector) // умножение матрицы на вектор 
                                                         // при диагональном предобуславливании,
                                                         // сразу применен ход для решения СЛАУ
        {
            Vector product = new Vector(n);

            for (uint i = 0; i < n; i++)
                product.vec[i] = 1 / Math.Sqrt(di.vec[i]) * vector.vec[i];

            return product;
        }

        public void LOS() // Локально-оптимальная схема без предобуславливания
        {
            uint index;
            real alpha, beta;

            Vector tmp = new Vector(n);

            r = pr - Mult(x);

            Array.Copy(r.vec, z.vec, n);
            p = Mult(z);

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

            Vector tmp = new Vector(n);

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

            Vector tmp = new Vector(n);

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

        public void CGM() // Метод сопряженных градиентов
        {
            uint index;
            real alpha, beta;
            real tmp;

            Vector temp = new Vector(n);

            r = pr - Mult(x);
            Array.Copy(r.vec, z.vec, n);

            for (index = 0; index < maxIter && (squareNorm = CalcNorm(r) / normPr) >= eps; index++)
            {
                temp = Mult(z);
                alpha = (r * r) / (temp * z);
                x = x + alpha * z;
                tmp = r * r;
                r = r - alpha * temp;
                beta = (r * r) / tmp;
                z = r + beta * z;
            }

            countIter = index;
        }

        public void CGMWithDi() // Метод сопряженных градиентов с диагональным предобуславливанием
        {
            uint index;
            real alpha, beta;
            real tmp;

            Vector fstTemp = new Vector(n);
            Vector sndTemp = new Vector(n);

            r = pr - Mult(x);
            z = MultDi(r);

            for (index = 0; index < maxIter && (squareNorm = CalcNorm(r) / normPr) >= eps; index++)
            {
                tmp = MultDi(r) * r;
                sndTemp = Mult(z);
                alpha = tmp / (sndTemp * z);
                x = x + alpha * z;
                r = r - alpha * sndTemp;
                fstTemp = MultDi(r);
                beta = (fstTemp * r) / tmp;
                z = fstTemp + beta * z;
            }

            countIter = index;
        }

        public void CGMWithCholesky() // Метод сопряженных градиентов с неполной факторизацией
        {
            uint index;
            real alpha, beta;
            real tmp;

            Vector fstTemp = new Vector(n);
            Vector sndTemp = new Vector(n);

            Cholesky();

            r = pr - Mult(x);
            z = MoveForCholesky(r);

            for (index = 0; index < maxIter && (squareNorm = CalcNorm(r) / normPr) >= eps; index++)
            {
                tmp = MoveForCholesky(r) * r;
                sndTemp = Mult(z);
                alpha = tmp / (sndTemp * z);
                x = x + alpha * z;
                r = r - alpha * sndTemp;
                fstTemp = MoveForCholesky(r);
                beta = (fstTemp * r) / tmp;
                z = fstTemp + beta * z;
            }

            countIter = index;
        }

        public void CGMAssymetric() // Метод сопряженных градиентов для несимметричных матриц
        {
            uint index;
            real alpha, beta;
            real tmp;

            Vector temp = new Vector(n);

            LU();

            r = DirectT(MultT(ReverseT(Direct(pr - Mult(x)))));
            Array.Copy(r.vec, z.vec, n);

            for (index = 0; index < maxIter && (squareNorm = CalcNorm(r) / normPr) >= eps; index++)
            {
                tmp = r * r;
                temp = DirectT(MultT(ReverseT(Direct(Mult(Reverse(z))))));
                alpha = tmp / (temp * z);
                x = x + alpha * z;
                r = r - alpha * temp;
                beta = (r * r) / tmp;
                z = r + beta * z;
            }

            countIter = index;

            x = Reverse(x);
        }

        public void CGMAssymetricDi() // Метод сопряженных градиентов для несимметричных матриц с диагональным предобуславливанием
        {
            uint index;
            real alpha, beta;
            real tmp;

            Vector temp = new Vector(n);

            r = MultDi(MultT(MultDi(MultDi(pr - Mult(x)))));

            Array.Copy(r.vec, z.vec, n);

            for (index = 0; index < maxIter && (squareNorm = CalcNorm(r) / normPr) >= eps; index++)
            {
                tmp = r * r;
                temp = MultDi(MultT(MultDi(MultDi(Mult(MultDi(z))))));
                alpha = tmp / (temp * z);
                x = x + alpha * z;
                r = r - alpha * temp;
                beta = (r * r) / tmp;
                z = r + beta * z;
            }

            countIter = index;

            x = MultDi(x);
        }

        private real CalcNorm(Vector vector)
        {
            real result = 0;

            for (uint i = 0; i < vector.vec.Length; i++)
                result += vector.vec[i] * vector.vec[i];

            return Math.Sqrt(result);
        }

        private void Cholesky()
        {
            real suml = 0;
            real sumdi = 0;

            for (uint i = 0; i < n; i++)
            {
                uint i0 = ig[i];
                uint i1 = ig[i + 1];

                for (uint k = i0; k < i1; k++)
                {
                    uint j = jg[k];
                    uint j0 = ig[j];
                    uint j1 = ig[j + 1];
                    uint ik = i0;
                    uint kj = j0;

                    while (ik < k && kj < j1)
                    {
                        if (jg[ik] == jg[kj])
                        {
                            suml += gglnew.vec[ik] * gglnew.vec[kj];
                            ik++;
                            kj++;
                        }
                        else
                        {
                            if (jg[ik] > jg[kj])
                                kj++;
                            else
                                ik++;
                        }
                    }

                    gglnew.vec[k] = (gglnew.vec[k] - suml) / dinew.vec[j];
                    sumdi += gglnew.vec[k] * gglnew.vec[k];
                    suml = 0;
                }

                dinew.vec[i] = Math.Sqrt(dinew.vec[i] - sumdi);
                sumdi = 0;
            }
        }

        private void LU()
        {
            real suml = 0;
            real sumu = 0;
            real sumdi = 0;

            for (uint i = 0; i < n; i++)
            {
                uint i0 = ig[i];
                uint i1 = ig[i + 1];

                for (uint k = i0; k < i1; k++)
                {
                    uint j = jg[k];
                    uint j0 = ig[j];
                    uint j1 = ig[j + 1];
                    uint ik = i0;
                    uint kj = j0;

                    while (ik < k && kj < j1)
                    {
                        if (jg[ik] == jg[kj])
                        {
                            suml += gglnew.vec[ik] * ggunew.vec[kj];
                            sumu += ggunew.vec[ik] * gglnew.vec[kj];
                            ik++;
                            kj++;
                        }
                        else
                        {
                            if (jg[ik] > jg[kj])
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

        private Vector DirectT(Vector vector) // U^-T
        {
            Vector result = new Vector(n);

            Array.Copy(vector.vec, result.vec, n);

            real sum = 0;

            for (uint i = 0; i < n; i++)
            {
                uint i0 = ig[i];
                uint i1 = ig[i + 1];

                for (uint k = i0; k < i1; k++)
                    sum += result.vec[jg[k]] * ggunew.vec[k];

                result.vec[i] -= sum;
                sum = 0;
            }

            return result;
        }

        private Vector ReverseT(Vector vector) // L^-T
        {
            Vector result = new Vector(n);

            Array.Copy(vector.vec, result.vec, n);

            for (int i = (int)n - 1; i >= 0; i--)
            {
                uint i0 = ig[i];
                uint i1 = ig[i + 1];
                result.vec[i] /= dinew.vec[i];

                for (uint k = i0; k < i1; k++)
                    result.vec[jg[k]] -= gglnew.vec[k] * result.vec[i];
            }

            return result;
        }

        private Vector MoveForCholesky(Vector vector)
        {

            Vector y = new Vector(n);
            Vector x = new Vector(n);

            Array.Copy(vector.vec, y.vec, n);

            real sum = 0;

            for (uint i = 0; i < n; i++) // Прямой ход
            {
                uint i0 = ig[i];
                uint i1 = ig[i + 1];

                for (uint k = i0; k < i1; k++)
                    sum += gglnew.vec[k] * y.vec[jg[k]];

                y.vec[i] = (y.vec[i] - sum) / dinew.vec[i];
                sum = 0;
            }

            Array.Copy(y.vec, x.vec, y.vec.Length);

            for (int i = (int)n - 1; i >= 0; i--) // Обратный ход
            {
                uint i0 = ig[i];
                uint i1 = ig[i + 1];
                x.vec[i] = y.vec[i] / dinew.vec[i];

                for (uint k = i0; k < i1; k++)
                    y.vec[jg[k]] -= gglnew.vec[k] * x.vec[i];
            }

            return x;
        }

        private Vector Direct(Vector vector)
        {
            Vector y = new Vector(n);

            Array.Copy(vector.vec, y.vec, n);

            real sum = 0;

            for (uint i = 0; i < n; i++)
            {
                uint i0 = ig[i];
                uint i1 = ig[i + 1];

                for (uint k = i0; k < i1; k++)
                    sum += gglnew.vec[k] * y.vec[jg[k]];

                y.vec[i] = (y.vec[i] - sum) / dinew.vec[i];
                sum = 0;
            }

            return y;
        }

        private Vector Reverse(Vector vector)
        {
            Vector result = new Vector(n);

            Array.Copy(vector.vec, result.vec, n);

            for (int i = (int)n - 1; i >= 0; i--)
            {
                uint i0 = ig[i];
                uint i1 = ig[i + 1];

                for (uint k = i0; k < i1; k++)
                    result.vec[jg[k]] -= ggunew.vec[k] * result.vec[i];
            }

            return result;
        }

        public void WriteToFile(string path, string str, real time)
        {
            using (var sw = new StreamWriter(path, true))
            {
                sw.WriteLine($"Method: {str}\n\n iterations: {countIter}, residual: {squareNorm.ToString("0.00E+0")}, time(ms): {time}\n");

                // for (uint i = 0; i < n; i++)
                //     sw.WriteLine(x.vec[i]);

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
            Array.Resize(ref ig, (int)n + 1);
            Array.Resize(ref jg, (int)(n * (n - 1) / 2));
            Array.Resize(ref di.vec, (int)n);
            Array.Resize(ref ggl.vec, (int)(n * (n - 1) / 2));
            Array.Resize(ref ggu.vec, (int)(n * (n - 1) / 2));
            Array.Resize(ref pr.vec, (int)n);
            Array.Resize(ref dinew.vec, (int)n);
            Array.Resize(ref gglnew.vec, (int)(n * (n - 1) / 2));
            Array.Resize(ref ggunew.vec, (int)(n * (n - 1) / 2));
            Array.Resize(ref r.vec, (int)n);
            Array.Resize(ref z.vec, (int)n);
            Array.Resize(ref p.vec, (int)n);
            Array.Resize(ref x.vec, (int)n);
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
            ig[0] = 0;

            Resize();

            for (uint i = 0; i < dim; i++)
            {
                di.vec[i] = 1.0 / (2 * (i + 1) - 1);
                ig[i + 1] = ig[i] + i;
            }

            for (uint i = 0; i < dim; i++)
            {
                real sum = 0;
                uint i0 = ig[i];
                uint i1 = ig[i + 1];
                uint j = i - (i1 - i0);

                for (uint k = i0; k < i1; k++, j++)
                {
                    sum = 1.0 / (i + 1 + j + 1 - 1);
                    ggl.vec[k] = sum;
                    ggu.vec[k] = sum;
                    jg[k] = j;
                }
            }

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