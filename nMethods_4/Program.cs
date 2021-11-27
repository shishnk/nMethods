using nMethods_4;

SNE sne = new("parametrs.txt", Test.sixthTest, Derivative.Analytic);
sne.MethodNewton();
sne.WriteToFile("x.txt");