using nMethods_6;

Spline spline = new("elements.txt", "points.txt", "parameters.txt");
spline.Compute();
spline.WriteToFile("result.txt");