namespace nMethods_6;

public record struct FiniteElement
{
    public double LeftBorder { get; init; }
    public double RightBorder { get; init; }
    public double Lenght { get; init; }

    public FiniteElement(double left, double right)
    {
        LeftBorder = left;
        RightBorder = right;
        Lenght = Math.Abs(RightBorder - LeftBorder);
    }

    public bool Contain(Point2D point)
    {
        if (point.X >= LeftBorder && point.X <= RightBorder)
            return true;
        else
            return false;
    }

    public static FiniteElement Parse(string elements)
    {
        var data = elements.Split(" ");
        FiniteElement element = new(double.Parse(data[0]), double.Parse(data[1]));

        return element;
    }

    public override string ToString()
    {
        return $"Element interval is [{LeftBorder}, {RightBorder}]";
    }
}