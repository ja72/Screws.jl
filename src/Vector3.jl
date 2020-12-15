
    """
    Cartesian 3×1 vector

        Vector3()
        Vector3(x,y,z)
        Vector3(vector)

    """
    const Vector3 = SVector{3, Float64}

    Vector3() = Vector3(fill(0.0, 3)) 

    export ô, î, ĵ, k̂
    """
    Origin vector ``\\hat{o}=\\pmatrix{0\\\\0\\\\0}``

        const ô = Vector3(0,0,0)

    """
    const ô = Vector3(0,0,0)
    """
    X-axis basis vector ``\\hat{\\imath}=\\pmatrix{1\\\\0\\\\0}``

        const î = Vector3(1,0,0)

    """
    const î = Vector3(1,0,0)
    """
    Y-axis basis vector ``\\hat{\\jmath}=\\pmatrix{0\\\\1\\\\0}``

        const ĵ = Vector3(0,1,0)

    """
    const ĵ = Vector3(0,1,0)
    """
    Z-axis basis vector ``\\hat{k}=\\pmatrix{0\\\\0\\\\1}``

        const k̂ = Vector3(0,0,1)

    """
    const k̂ = Vector3(0,0,1)

    function  LinearAlgebra.norm(a::Vector3, p::Real=2)
        return LinearAlgebra.norm(a.data, p)
    end

    """
    Returns a normalized a vector whose magnitude is 1.0 or 0.0

        n = normalize(v)

    """
    function LinearAlgebra.normalize(v::Vector3)::Vector3
        m2 = dot(v,v)
        if m2>0 && m2!=1.0
            return v/sqrt(m2)
        else
            return v
        end
    end

    """
    Dot product between vectors ``\\pmatrix{x \\\\ y \\\\ z} \\cdot \\pmatrix{u \\\\ v \\\\ w} = x*u + y*v + z*w``

        dot(a::Vector3, b::Vector3)
        a ⋅ b => dot(a,b)
        
    """
    @inline dot(a::Vector3, b::Vector3)::Float64 = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]

    """
    Cross product of vectors

        cross(a::Vector3, b::Vector3)::Vector3
        a × b => cross(a,b)

    """
    function cross(a::Vector3, b::Vector3)::Vector3 
        Vector3(a[2]*b[3]-a[3]*b[2],
            a[3]*b[1]-a[1]*b[3],
            a[1]*b[2]-a[2]*b[1])
    end

    """
    Cross product matrix operator ``\\pmatrix{x \\\\ y \\\\ z}\\times =  \\pmatrix{0 & -z & y \\\\ z & 0 & -x \\\\ -y & x & 0}``

        cross(a::Vector3)
        ×(a::Vector3)

    """
    cross(a::Vector3)::Matrix3 = Matrix3(0,-a[3],a[2],a[3],0,-a[1],-a[2],a[1],0)

    """
    Double cross product matrix operator ``\\pmatrix{x \\\\ y \\\\ z}\\times \\pmatrix{x \\\\ y \\\\ z}\\times =  \\pmatrix{-y^2-z^2 & x y & x z \\\\ x y & -x^2-z^2 & y z \\\\ x z & y z & -x^2-z^2}``

        cross2(a::Vector3)

    """
    function cross2(a::Vector3)::Matrix3 
        return Matrix3( 
            -a[2]^2-a[3]^2, a[1]*a[2], a[1]*a[3],  
            a[1]*a[2], -a[1]^2-a[3]^2, a[2]*a[3],
            a[1]*a[3], a[2]*a[3], -a[1]^2-a[2]^2)
    end

    function outer(a::Vector3, b::Vector3)::Matrix3
        return Matrix3(
            a[1]*b[1], a[1]*b[2], a[1]*b[3],
            a[2]*b[1], a[2]*b[2], a[2]*b[3],
            a[3]*b[1], a[3]*b[2], a[3]*b[3])
    end

    function Base.angle(a::Vector3, b::Vector3)::Float64
        x = dot(a,b)
        y = cross(a,b).Magnitude
        return atan(y, x)
    end

    """
    Rotational interpolation between two vectors
    
        slerp(a,b,t)    # where t=0..1

    """
    function slerp(a::Vector3, b::Vector3, t::Float64) :: Vector3
        θ = angle(a, b)
        y = sin(θ)
        return sin((1-t)*θ)/y*a + sin(t*θ)/y*b
    end

    function Base.getproperty(v::Vector3, name::Symbol)
        if name===:X  
            return v[1]
        end
        if name===:Y  
            return v[2]
        end
        if name===:Z
            return v[3]
        end
        if name===:SumSquares
            return dot(v,v)
        end
        if name===:Magnitude
            return sqrt(dot(v,v))
        end
        return getfield(v, name)
    end

    function Base.propertynames(v::Vector3)
        return (:X, :Y, :Z, :Magnitude, :SumSquares)
    end
