global using real = System.Double;

namespace nMethods_4
{
    public class SNE
    {
        private real eps; // точность решения СНУ
        private real maxIter; // максимальное кол-во итераций
        private real[] x; // начальное приближение
        private real[] f; // вектор правой части
        private real[,] matrixJacobi;

        public SNE(string path)
        {
            try
            { // считывание данных
                using (var sr = new StreamReader(path))
                {
                    eps = real.Parse(sr.ReadLine());
                    maxIter = real.Parse(sr.ReadLine());
                    x = sr.ReadLine().Split(" ").Select(value => real.Parse(value)).ToArray();
                }

                matrixJacobi = new real[4, 4];
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.Message);
            }
        }

        private void CalcElementsJacobi()
        {
            matrixJacobi[0, 0] = 1;
            matrixJacobi[0, 1] = 1;
            matrixJacobi[1, 0] = 2 * x[0];
            matrixJacobi[1, 1] = 2 * x[1];
        }

        private void CalcElementsF()
        {

        }

        private void MethodGauss()
        {
            real[] result = new real[matrixJacobi.GetLength(0)];
            real max;
            uint k, index;

            k = 0;

            while (k < matrixJacobi.GetLength(0))
            {
                max = Math.Abs(matrixJacobi[k,k]);
                index = k;

                for (uint i = 0; i < matrixJacobi.GetLength(0); i++)
                {
                    if (Math.Abs(matrixJacobi[i,k]) > max)
                    {
                        max = Math.Abs(matrixJacobi[i,k]);
                        index = i;
                    }
                }

                for (uint j = 0; j < matrixJacobi.GetLength(0); j++)
                {
                    real temp = matrixJacobi[k,j];
                    matrixJacobi[k, j] = matrixJacobi[index, j];
                    matrixJacobi[index, j] = temp;
                }

                real tmp = f[k];
                f[k] = f[index];
                f[index] = tmp; // TODO
            }

        }
    }
}
