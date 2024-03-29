﻿using System;
using System.Diagnostics;

namespace nMethods_3
{
    class Program
    {
        static void Main(string[] args)
        {
            SLAE slau = new SLAE("kuslau.txt", "ig.txt", "jg.txt", "ggl.txt",
                                 "ggu.txt", "di.txt", "pr.txt");

            var watch = new Stopwatch();

            #region LOS
            // watch.Start();
            // slau.LOS();
            // watch.Stop();

            // slau.WriteToFile("x.txt", "LOS", watch.ElapsedMilliseconds * 1000);
            // slau.Clear();

            // watch.Start();
            // slau.LOSWithLU();
            // watch.Stop();

            // slau.WriteToFile("x.txt", "LOS + LU", watch.ElapsedMilliseconds * 1000);
            // slau.Clear();

            // watch.Start();
            // slau.LOSWithDi();
            // watch.Stop();

            // slau.WriteToFile("x.txt", "LOS + DI", watch.ElapsedMilliseconds * 1000);
            #endregion

            #region Hilbert
            // Console.WriteLine("Введите размерность:");

            // slau.GenMatrixHilbert(uint.Parse(Console.ReadLine()));
            // watch.Start();
            // slau.LOS();
            // watch.Stop();

            // slau.WriteToFile("x(Hilbert).txt", "LOS | Hilbert", watch.ElapsedMilliseconds * 1000);
            // slau.Clear();

            // watch.Start();
            // slau.LOSWithLU();
            // watch.Stop();

            // slau.WriteToFile("x(Hilbert).txt", "LOS + LU | Hilbert", watch.ElapsedMilliseconds * 1000);
            // slau.Clear();

            // watch.Start();
            // slau.LOSWithDi();
            // watch.Stop();

            // slau.WriteToFile("x(Hilbert).txt", "LOS + DI | Hilbert", watch.ElapsedMilliseconds * 1000);
            #endregion

            #region CMG
            // watch.Start();
            // slau.CGM();
            // watch.Stop();

            // slau.WriteToFile("x.txt", "CMG", watch.ElapsedMilliseconds * 1000);
            // slau.Clear();

            // watch.Start();
            // slau.CGMWithCholesky();
            // watch.Stop();

            // slau.WriteToFile("x.txt", "CMG + LL^T", watch.ElapsedMilliseconds * 1000);
            // slau.Clear();

            // watch.Start();
            // slau.CGMWithDi();
            // watch.Stop();

            // slau.WriteToFile("x.txt", "CMG + DI", watch.ElapsedMilliseconds * 1000);
            // slau.Clear();            
            #endregion

            #region CMG for asymmetric matrixes
            watch.Start();
            slau.CGMAssymetric();
            watch.Stop();

            slau.WriteToFile("x.txt", "CMG for assymetric matrixes", watch.ElapsedMilliseconds * 1000);
            slau.Clear();

            watch.Start();
            slau.CGMAssymetricDi();
            watch.Stop();

            slau.WriteToFile("x.txt", "CMG for assymetric matrixes + DI", watch.ElapsedMilliseconds * 1000);
            slau.Clear();
            #endregion
        }
    }
}