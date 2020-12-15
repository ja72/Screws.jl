# Screws.jl

Definitions for Twists and Wrenches for use in 3D rigid body dynamics.

## Types/Aliases

 - Vector3 = SVector{3, Float64}
 - Matrix3 = SMatrix{3, 3, Float64}
 - Matrix6 = SMatrix{2,2, Matrix3}
 - Twist(value,position,pitch)
 - Wrench(value,position,pitch)
 - SpatialInertia(ll,la,al,aa)
 - SpatialMobility(ll,la,al,aa)

## Functions/Methods

 - spi(m,Ic,cg)
 - spm(m,Ic,cg)
 - dot(v,f), dot(f,v)
