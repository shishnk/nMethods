using nMethods_4;

SNE sne = new("parametrs.txt", Test.firstTest);
sne.MethodNewton();
sne.WriteToFile("x.txt");