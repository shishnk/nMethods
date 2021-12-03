namespace nMethods_4
{
    public class Vector
    {
        public real[] vec;
        public int Length { get; init; }

        public Vector(int dim)
        {
            vec = new real[dim];
            Length = dim;
        }

        public real this[int index]
        {
            get => vec[index];
            set => vec[index] = value;
        }

        public static Vector operator *(real constant, Vector vector)
        {
            Vector result = new(vector.Length);

            for (int i = 0; i < vector.Length; i++)
                result[i] = constant * vector[i];

            return result;
        }

        public static Vector operator +(Vector fstVector, Vector sndVector)
        {
            Vector result = new(fstVector.Length);

            for (int i = 0; i < fstVector.Length; i++)
            {
                result[i] = fstVector[i] + sndVector[i];
            }

            return result;
        }

        public static Vector operator -(Vector vector)
        {
            Vector result = new(vector.Length);

            for (int i = 0; i < vector.Length; i++)
                result[i] = vector[i] * (-1);

            return result;
        }

        public static void Copy(Vector source, Vector destination)
        {
            for (int i = 0; i < source.Length; i++)
                destination[i] = source[i];
        }

        public static void Clear(Vector vector)
        {
            for (int i = 0; i < vector.Length; i++)
                vector[i] = 0;
        }

        public static void Resize(ref Vector vector, int lenght)
        {
            vector = new(lenght);
        }
    }
}