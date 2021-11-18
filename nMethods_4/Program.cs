namespace nMethods_4
{
    class Program
    {
        public static void Main(string[] args)
        {
            SNE sne = new("parametrs.txt");
            sne.MethodNewton();
            sne.WriteToFile("x.txt");
        }
    }
}