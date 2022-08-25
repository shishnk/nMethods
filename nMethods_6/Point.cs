namespace nMethods_6;

public readonly record struct Point(double X, double Y)
{
    public static Point operator +(Point point, (double, double) value)
         => new(point.X + value.Item1, point.Y + value.Item2);


    public override string ToString()
    {
        return $"({X},{Y})";
    }
}