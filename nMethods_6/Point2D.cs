namespace nMethods_6;

public record struct Point2D
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

    public static Point2D operator +(Point2D point, (double, double) value)
         => new(point.X + value.Item1, point.Y + value.Item2);

    //public Point2D ProjectionX() => new(X, 0);

    //public Point2D ProjectionY() => new(0, Y);

    public override string ToString()
    {
        return $"({X},{Y})";
    }
}