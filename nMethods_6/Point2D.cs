namespace nMethods_6;
public struct Point2D
{
    public double X { get; init; }
    public double Y { get; init; }

    public Point2D(double x, double y)
    {
        X = x;
        Y = y;
    }

    public static Point2D Parse(string points)
    {
        var data = points.Split(" ");
        Point2D point = new(double.Parse(data[0]), double.Parse(data[1]));

        return point;
    }

    public override string ToString()
    {
        return $"({X},{Y})";
    }
}