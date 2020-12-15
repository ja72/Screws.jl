
    """
    Cartesian 3×3 matrix

        Matrix3()
        Matrix3(a11,a12,a13,a21,a23,a23,a31,a32,a33)*
        Matrix3(matrix)

    (*) definition is row major, overriding the column major definition of `SMatrix{}`

    """
    const Matrix3 = SMatrix{3, 3, Float64}

    Matrix3() = Matrix3(fill(0.0, 3,3))
    # Define matrix elements in row major order
    Matrix3(a11,a12,a13,a21,a22,a23,a31,a32,a33) = Matrix3([a11 a12 a13; a21 a22 a23; a31 a32 a33])

    export Ô, Î

    """
    Zero 3×3 matrix ``\\hat{O} = \\pmatrix{0 & 0 & 0 \\\\ 0 & 0 & 0 \\\\ 0 & 0 & 0}``

        const Ô = Matrix3(
            0,0,0,
            0,0,0,
            0,0,0)

    """
    const Ô = Matrix3()
    """
    Identity 3×3 matrix ``\\hat{I} = \\pmatrix{1 & 0 & 0 \\\\ 0 & 1 & 0 \\\\ 0 & 0 & 1}``

        const Î = Matrix3(
            1,0,0,
            0,1,0,
            0,0,1)

    """
    const Î = Matrix3(I)

    function Base.inv(m::Matrix3)::Matrix3
        t2 = m[1,1]*m[2,2]*m[3,3]
        t3 = m[1,2]*m[2,3]*m[3,1]
        t4 = m[1,3]*m[2,1]*m[3,2]
        t7 = m[1,1]*m[2,3]*m[3,2]
        t8 = m[1,2]*m[2,1]*m[3,3]
        t9 = m[1,3]*m[2,2]*m[3,1]
        t6 = 1.0/(t2+t3+t4-t7-t8-t9)
        return Matrix3(
             t6*(m[2,2]*m[3,3]-m[2,3]*m[3,2]),
            -t6*(m[1,2]*m[3,3]-m[1,3]*m[3,2]),
             t6*(m[1,2]*m[2,3]-m[1,3]*m[2,2]),
            -t6*(m[2,1]*m[3,3]-m[2,3]*m[3,1]),
             t6*(m[1,1]*m[3,3]-m[1,3]*m[3,1]),
            -t6*(m[1,1]*m[2,3]-m[1,3]*m[2,1]),
             t6*(m[2,1]*m[3,2]-m[2,2]*m[3,1]),
            -t6*(m[1,1]*m[3,2]-m[1,2]*m[3,1]),
             t6*(m[1,1]*m[2,2]-m[1,2]*m[2,1]))          
    end

    function solve(m::Matrix3, v::Vector3)::Vector3
        t = (-m[1,2]*m[3,3]+m[1,3]*m[3,2])*m[2,1]+(-m[2,2]*m[1,3]+m[2,3]*m[1,2])*m[3,1]+(m[2,2]*m[3,3]-m[2,3]*m[3,2])*m[1,1]
        return Vector3(
            (v[1]*( m[2,2]*m[3,3] - m[2,3]*m[3,2]) + v[2]*(-m[1,2]*m[3,3] + m[1,3]*m[3,2]) + v[3]*(-m[2,2]*m[1,3] + m[2,3]*m[1,2]))/t,
            (v[1]*(-m[2,1]*m[3,3] + m[2,3]*m[3,1]) + v[2]*( m[1,1]*m[3,3] - m[1,3]*m[3,1]) + v[3]*( m[2,1]*m[1,3] - m[2,3]*m[1,1]))/t,
            (v[1]*( m[2,1]*m[3,2] - m[2,2]*m[3,1]) + v[2]*(-m[1,1]*m[3,2] + m[1,2]*m[3,1]) + v[3]*(-m[2,1]*m[1,2] + m[2,2]*m[1,1]))/t) 
    end

    function Base.:\(A::Matrix3, b::Vector3)::Vector3
        return solve(A,b)
    end    

    function Base.getproperty(m::Matrix3, name::Symbol)
        if name===:A11; return m[1,1]; end
        if name===:A12; return m[1,2]; end
        if name===:A13; return m[1,3]; end
        if name===:A21; return m[2,1]; end
        if name===:A22; return m[2,2]; end
        if name===:A23; return m[2,3]; end
        if name===:A31; return m[3,1]; end
        if name===:A32; return m[3,2]; end
        if name===:A33; return m[3,3]; end
        if name===:Trace; return m[1,1]+m[2,2]+m[3,3]; end
        if name===:Determinat; return m[1,1]*(m[2,2]*m[3,3]-m[2,3]*m[3,2])+m[1,2]*(m[2,3]*m[3,1]-m[2,1]*m[3,3])+m[1,3]*(m[2,1]*m[3,2]-m[2,2]*m[3,1]); end
        if name===:Row1; return Vector3(m[1,1], m[1,2], m[1,3]); end
        if name===:Row2; return Vector3(m[2,1], m[2,2], m[2,3]); end
        if name===:Row3; return Vector3(m[3,1], m[3,2], m[3,3]); end
        if name===:Column1; return Vector3(m[1,1], m[2,1], m[3,1]); end
        if name===:Column2; return Vector3(m[1,2], m[2,2], m[3,2]); end
        if name===:Column3; return Vector3(m[1,3], m[2,3], m[3,3]); end
        return getfield(m, name)
    end

    function Base.propertynames(m::Matrix3)
        return (:A11, :A12, :A13, :A21, :A22, :A23, :A31, :A32, :A33,
        :Row1, :Row2, :Row3, :Column1, :Column2, :Column3, :Trace, :Determinat)
    end
