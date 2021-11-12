using real = System.Double;

namespace nMethods_3
{
    class Vector 
    {
      public real[] vec;

    public Vector(uint dim)
    {
        vec = new real[dim];
    }

    public static real operator *(Vector fstVec, Vector scdVec)
    {
        real result = 0;

        for (uint i = 0; i < fstVec.vec.Length; i++)
            result += fstVec.vec[i] * scdVec.vec[i];

        return result;
    }

    public static Vector operator *(real constant, Vector vector)
    {
        Vector result = new Vector((uint)vector.vec.Length);

        for (uint i = 0; i < vector.vec.Length; i++)
            result.vec[i] = vector.vec[i] * constant;

        return result;
    }

    public static Vector operator +(Vector fstVec, Vector scdVec)
    {
        Vector result = new Vector((uint)fstVec.vec.Length);

        for (uint i = 0; i < fstVec.vec.Length; i++)
            result.vec[i] = fstVec.vec[i] + scdVec.vec[i];

        return result;
    }

    public static Vector operator -(Vector fstVec, Vector scdVec)
    {
        Vector result = new Vector((uint)fstVec.vec.Length);

        for (uint i = 0; i < fstVec.vec.Length; i++)
            result.vec[i] = fstVec.vec[i] - scdVec.vec[i];

        return result;
    }
}
}