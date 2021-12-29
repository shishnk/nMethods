namespace nMethods_6;
public struct FiniteElement
{
    public double leftBorder { get; init; }
    public double rightBorder { get; init; }
    public double Lenght { get; init; }

    public FiniteElement(double left, double right)
    {
        leftBorder = left;
        rightBorder = right;
        Lenght = Math.Abs(rightBorder - leftBorder);
    }

    public bool Contain(Point2D point)
    {
        if (point.X >= leftBorder && point.X <= rightBorder)
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
        return $"Element interval is [{leftBorder}, {rightBorder}]";
    }
}