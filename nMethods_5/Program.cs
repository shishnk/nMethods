using nMethods_5;

Eigenvalue eigenvalues = new(new Matrix("matrix.txt"));
eigenvalues.Compute();
eigenvalues.WriteToFile("result.txt");