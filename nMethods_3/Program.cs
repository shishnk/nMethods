using System.Diagnostics;

namespace nMethods_3
{
    class Program
    {
        static void Main(string[] args)
        {
            SLAU slau = new SLAU("kuslau.txt", "ig.txt", "jg.txt", "ggl.txt",
                                 "ggu.txt", "di.txt", "pr.txt");

            var watch = new Stopwatch();

            watch.Start();
            slau.LOS();
            watch.Stop();
            
            slau.WriteToFile("x.txt", "LOS",watch.ElapsedMilliseconds * 1000);
            slau.Clear();

            watch.Start();
            slau.LOSWithLU();
            watch.Stop();

            slau.WriteToFile("x.txt", "LOS + LU",watch.ElapsedMilliseconds * 1000);
            slau.Clear();

            watch.Start();
            slau.LOSWithDi();
            watch.Stop();

            slau.WriteToFile("x.txt", "LOS + DI", watch.ElapsedMilliseconds * 1000);
        }
    }
}