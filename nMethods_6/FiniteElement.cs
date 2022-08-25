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

    public bool Contain(Point point)
        => (point.X >= LeftBorder && point.X <= RightBorder) ? true : false;

    public override string ToString()
    {
        return $"Element interval is [{LeftBorder}, {RightBorder}]";
    }
}