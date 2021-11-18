namespace nMethods_4
{
    public class Vector
    {
        public real[] vec;

        public Vector(int dim)
        {
            vec = new real[dim];
        }

        public static Vector operator *(real constant, Vector vector)
        {
            Vector result = new(vector.vec.Length);

            for (uint i = 0; i < vector.vec.Length; i++)
                result.vec[i] = constant * vector.vec[i];

            return result;
        }

        public static Vector operator +(Vector fstVector, Vector sndVector)
        {
            Vector result = new(fstVector.vec.Length);

            for (uint i = 0; i < fstVector.vec.Length; i++)
            {
                result.vec[i] = fstVector.vec[i] + sndVector.vec[i];
            }

            return result;
        }

        public static Vector operator -(Vector vector)
        {
            Vector result = new(vector.vec.Length);

            for (uint i = 0; i < vector.vec.Length; i++)
                result.vec[i] = vector.vec[i] * (-1);
            
            return result;
        }
    }
}