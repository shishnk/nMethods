using real = System.Double;

namespace nMethods_3
{
    class Vector<T>
    {
      public T[] vec;

        public Vector(uint dim)
        {
            vec = new T[dim];
        }

        public static real operator *(Vector<T> fstVec, Vector<T> scdVec)
        {
            real result = 0;

            for (uint i = 0; i < fstVec.vec.Length; i++)
                result += (real)(object)fstVec.vec[i] * (real)(object)scdVec.vec[i];

            return result;
        }

        public static Vector<real> operator *(real constant, Vector<T> vector)
        {
            Vector<real> result = new Vector<real>((uint)vector.vec.Length);

            for (uint i = 0; i < vector.vec.Length; i++)
                result.vec[i] = (real)(object)vector.vec[i] * constant;

            return result;
        }

        public static Vector<real> operator +(Vector<T> fstVec, Vector<T> scdVec)
        {
            Vector<real> result = new Vector<real>((uint)fstVec.vec.Length);

            for (uint i = 0; i < fstVec.vec.Length; i++)
                result.vec[i] = (real)(object)fstVec.vec[i] + (real)(object)scdVec.vec[i];

            return result;
        }

        public static Vector<real> operator -(Vector<T> fstVec, Vector<T> scdVec)
        {
            Vector<real> result = new Vector<real>((uint)fstVec.vec.Length);

            for (uint i = 0; i < fstVec.vec.Length; i++)
                result.vec[i] = (real)(object)fstVec.vec[i] - (real)(object)scdVec.vec[i];

            return result;
        }   
    }
}