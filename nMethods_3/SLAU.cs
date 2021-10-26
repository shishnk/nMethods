using System;
using System.IO;
using System.Linq;

using real = System.Double;

namespace nMethods_3
{
    class SLAU
    {
        private uint n;
        private uint maxIter;
        private uint countIter;
        private real eps;
        private real squareNorm;
        private Vector<uint> ig;
        private Vector<uint> jg;
        private Vector<real> ggl;
        private Vector<real> ggu;
        private Vector<real> di;
        private Vector<real> gglnew;
        private Vector<real> ggunew;
        private Vector<real> dinew;
        private Vector<real> pr;
        private Vector<real> r;
        private Vector<real> z;
        private Vector<real> p;
        private Vector<real> x;

        public SLAU(string pathParametrs, string pathRows, string pathColumns, string pathLower,
                    string pathUpper, string pathDiag, string pathVector)
        {
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

                Array.Copy(di.vec, dinew.vec, n);
                Array.Copy(ggl.vec, gglnew.vec, (uint)ggl.vec.Length);
                Array.Copy(ggu.vec, ggunew.vec, (uint)ggu.vec.Length);

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

        private Vector<real> MultDi(Vector<real> vector)
        {
            Vector<real> product = new Vector<real>(n);

            for (uint i = 0; i < n; i++)
                product.vec[i] = 1 / Math.Sqrt(di.vec[i]) * vector.vec[i];

            return product;
        }

        public void LOS()
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

        public void LOSWithLU()
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

        public void LOSWithDi() // TODO
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

            for (int i = (int)(n - 1); i >= 0; i--)
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
    }
}