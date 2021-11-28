using nMethods_4;

SNE sne = new("parametrs.txt", Test.fourthTest, Derivative.Numerical, ToSquareMatrix.Convolution);
sne.MethodNewton();
sne.WriteToFile("x.txt");