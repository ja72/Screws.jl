"""

Definitions for Twists and Wrenches for use in 3D rigid body dynamics

    $(EXPORTS)

    $(SIGNATURES)

    $(METHODLIST)

"""
module Screws
using LinearAlgebra
using StaticArrays
using Rotations
using DocStringExtensions
using Cartesian

using Base: @propagate_inbounds

import Base.convert, Base.*
import LinearAlgebra.cross, LinearAlgebra.dot
import Cartesian.outer

export 
    Twist,
    Wrench,
    Matrix6,
    SpatialInertia,
    SpatialMobility,
    spi,
    spm,
    dot

const Matrix6 = SMatrix{2,2, Matrix3}

for SpatialVector in (:Twist, :Wrench)
    @eval begin
        # Data Structures
        struct $SpatialVector
            linear::Vector3
            angular::Vector3
        
            @inline function $SpatialVector(linear::AbstractVector, angular::AbstractVector) where T
                @boundscheck( size(linear,1)==3 || throw(DimensionMismatch()))
                @boundscheck( size(angular,1)==3 || throw(DimensionMismatch()))
                return new(linear, angular)
            end
        end
        
        # Basics
        Base.eltype(::Type{$SpatialVector})  = Float64
        """
        $(SIGNATURES)
        Linear (first) part of the screw
        """        
        linear(v::$SpatialVector) = v.linear
        """
        $(SIGNATURES)
        Angular (second) part of the screw
        """        
        angular(v::$SpatialVector) = v.angular

        # Construct/convert given another spatial vector
        $SpatialVector(v::$SpatialVector) = v
        Base.convert(::Type{V}, v::$SpatialVector) where {V<:$SpatialVector} = V(v)

        # Construct/convert to SVector
        StaticArrays.SVector{6}(v::$SpatialVector) = convert(SVector{6}, [linear(v); angular(v)])
        StaticArrays.SVector(v::$SpatialVector) = SVector{6}(v)
        StaticArrays.SArray(v::$SpatialVector) = SVector(v)
        Base.convert(::Type{SA}, v::$SpatialVector) where {SA <: SArray} = SA(v)

        function Base.Array(v::$SpatialVector)
            Base.depwarn("Array(v) has been deprecated. Please use SVector(v) instead.", :Array)
            SVector(v)
        end
    end
end


"""
$(TYPEDEF)

A Twist represents a motion screw. A screw groups together a linear and an angular vector
into a single 6-vector. For twists the linear part transforms with location, and the angular
part remains constant.

"""
Twist

"""
$(TYPEDEF)

A Wrench represents a force screw. A screw groups together a linear and an angular vector
into a single 6-vector. For wrenches  the linear remains constant, while  the angular
part transforms with location.

"""
Wrench

"""
$(SIGNATURES)

Construct a [`Twist`](@ref) from a value vector, a position vector and a pitch

```math
\\mathrm{Twist}(\\boldsymbol{v}, \\boldsymbol{r}, h) = 
\\pmatrix{ \\boldsymbol{r}\\times \\boldsymbol{v} + h\\,\\boldsymbol{v} \\\\ \\boldsymbol{v} }
```
"""
Twist(value::Vector3, position::Vector3, pitch::Float64 ) = Twist( cross(position, value) + pitch*value, value)
#Twist(value::Vector3, position::Vector3) = Twist(value, position, 0.0)
"""
$(SIGNATURES)

Construct a pure Twist
```math
\\mathrm{Twist}(\\boldsymbol{v}) = 
\\pmatrix{ \\boldsymbol{v} \\\\ \\boldsymbol{0} }
```
"""
Twist(pure::Vector3 ) = Twist(pure, ô)
Twist() = Twist(ô)

function Base.:/(v::Twist, g::AbstractFloat)::Twist
    Wrench(v.linear/g, v.angular/g)
end


"""
$(SIGNATURES)
Construct a [`Wrench`](@ref) from a value vector, a position vector and a pitch

```math
\\mathrm{Wrench}(\\boldsymbol{v}, \\boldsymbol{r}, h) = 
\\pmatrix{ \\boldsymbol{v} \\\\ \\boldsymbol{r}\\times \\boldsymbol{v} + h\\,\\boldsymbol{v} }
```
"""
Wrench(value::Vector3, position::Vector3, pitch::Float64 ) = Wrench(value, cross(position, value) + pitch*value)
#Wrench(value::Vector3, position::Vector3) = Wrench(value, position, 0.0)

"""
$(SIGNATURES)

Construct a pure Wrench

```math
\\mathrm{Wrench}(\\boldsymbol{\tau}) = 
\\pmatrix{ \\boldsymbol{0} \\\\ \\boldsymbol{\tau}  }
```
"""
Wrench(pure::Vector3) = Wrench(ô, pure)
Wrench() = Wrench(ô)

value(v::Twist) = v.angular
axis(v::Twist) = v.linear
value(v::Wrench) = v.linear
axis(v::Wrench) = v.angular

function Base.:/(f::Wrench, g::AbstractFloat)::Wrench
    Wrench(f.linear/g, f.angular/g)
end

"""
$(SIGNATURES)

Dot products between Twists and Wrenches, used in power calculations.
"""
dot(v::Twist, f::Wrench) = dot(v.linear, f.linear) + dot(v.angular, f.angular)
dot(f::Wrench, v::Twist) = dot(v,f)

function outer(a::Twist, b::Twist)::Matrix6
    Matrix6(a.linear .* b.linear', a.linear .* b.angular',
    a.angular .* b.linear', a.angular .* b.angular')
end
function outer(a::Wrench, b::Wrench)::Matrix6
    Matrix6(a.linear .* b.linear', a.linear .* b.angular',
    a.angular .* b.linear', a.angular .* b.angular')
end


"""
$(SIGNATURES)

Cross products between Twist and Wrencehs, used in derivatives of rotating frames.
"""
cross(v::Twist, w::Twist) = Twist( cross(v.angular, w.linear) + cross(v.linear, w.angular), cross(v.angular, w.angular) )
cross(v::Twist, f::Wrench) = Wrench( cross(v.angular, w.linear), cross(v.linear, w.linear) + cross(v.angular, w.angular) )
cross(f::Wrench, v::Twist) = Twist( cross(f.linear, v.linear) + cross(f.angular, v.angular) + cross(f.linear, v.angular) )
cross(f::Wrench, g::Wrench) = Wrench( cross(f.linear, g.linear), cross(f.angular, g.linear) + cross(f.linear, g.angular) )

for SpatialMatrix in (:SpatialInertia, :SpatialMobility)
    @eval begin
        # Data Structures
        struct $SpatialMatrix
            ll::Matrix3
            la::Matrix3
            al::Matrix3
            aa::Matrix3

            @inline function $SpatialMatrix(ll::AbstractMatrix, la::AbstractMatrix, al::AbstractMatrix, aa::AbstractMatrix) where {T}
                @boundscheck( (size(ll,1)==3 && size(ll,2)==3) || throw(DimensionMismatch("$SpatialMatrix needs $ll to be a 3×3 matrix.")))
                @boundscheck( (size(la,1)==3 && size(la,2)==3) || throw(DimensionMismatch("$SpatialMatrix needs $la to be a 3×3 matrix.")))
                @boundscheck( (size(al,1)==3 && size(al,2)==3) || throw(DimensionMismatch("$SpatialMatrix needs $al to be a 3×3 matrix.")))
                @boundscheck( (size(aa,1)==3 && size(aa,2)==3) || throw(DimensionMismatch("$SpatialMatrix needs $aa to be a 3×3 matrix.")))
                # if !(size(ll,1) == size(ll,2) == 3) || !(size(la,1) == size(la,2) == 3) || !(size(al,1) == size(al,2) == 3) || !(size(aa,1) == size(aa,2) == 3)
                #     throw(DimensionMismatch("$SpatialMatrix needs four matrices of size 3×3"))
                # end
                return new(ll,la,al,aa)
            end
        end
        # Basics
        Base.eltype(::Type{$SpatialMatrix})  = Float64
        """
        $(SIGNATURES)
        LinearLinear (top-left) part of the matrix
        """        
        linearlinear(m::$SpatialMatrix) = m.ll
        """
        $(SIGNATURES)
        LinearAngular (top-right) part of the matrix
        """        
        linearangular(m::$SpatialMatrix) = m.la
        """
        $(SIGNATURES)
        AngularLinear (bottome-left) part of the matrix
        """        
        angularlinear(m::$SpatialMatrix) = m.al
        """
        $(SIGNATURES)
        AngularAngular (bottome-right) part of the matrix
        """        
        angularangular(m::$SpatialMatrix) = m.aa

        # Construct/convert given another spatial vector
        $SpatialMatrix(m::$SpatialMatrix) = m
        function $SpatialMatrix(m::Matrix6)
            $SpatialMatrix(m[1,1], m[1,2], m[2,1], m[2,2])
        end
        Base.convert(::Type{MType}, m::$SpatialMatrix) where {MType<:$SpatialMatrix} = MType(m)

        # Construct/convert to SVector
        StaticArrays.SMatrix{6,6}(m::$SpatialMatrix) = convert(SMatrix{6,6}, 
            [m.ll m.la; m.al m.aa])        
        StaticArrays.SMatrix(m::$SpatialMatrix) = SMatrix{6,6}(m)
        StaticArrays.SArray(m::$SpatialMatrix) = SMatrix(m)
        Base.convert(::Type{SA}, m::$SpatialMatrix) where {SA <: SArray} = SA(m)

        function Base.Array(m::$SpatialMatrix)
            Base.depwarn("Array(m) has been deprecated. Please use SMatrix(m) instead.", :Array)
            SMatrix(m)
        end
    end #eval
end #for


"""
$(TYPEDEF)


"""
SpatialInertia

"""
$(TYPEDEF)


"""
SpatialMobility


spi(m::AbstractFloat, Ic::Matrix3, cg::Vector3) = SpatialInertia( m*Î, -m*cross(cg), m*cross(cg), Ic-m*cross2(cg))
spm(m::AbstractFloat, Ic::Matrix3, cg::Vector3) = SpatialMobility( (1/m)*Î-cross(cg)*inv(Ic)*cross(cg), cross(cg)*inv(Ic), -inv(Ic)*cross(cg), inv(Ic))

function Base.:*(I::SpatialInertia, v::Twist)::Wrench
    Wrench(I.ll*v.linear + I.la*v.angular, I.al*v.linear + I.aa*v.angular)
end

function Base.:*(M::SpatialMobility, f::Wrench)::Twist
    Twist(M.ll*f.linear + M.la*f.angular, M.al*f.linear + M.aa*f.angular)    
end
    
end # module