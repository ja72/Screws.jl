# Source for Cartesian package. License is MIT https://github.com/ja72/Cartesian.jl/blob/master/LICENSE

"""

Definitions for cartesian vectors and matrices used in 3D mechanics. 
The repository is located at https://github.com/ja72/Cartesian.jl

    $(EXPORTS)

    $(SIGNATURES)

    $(METHODLIST)

"""
module Cartesian
using LinearAlgebra
using StaticArrays
using Rotations
using DocStringExtensions

using Base: @propagate_inbounds

import Base.convert, Base.*
import LinearAlgebra.cross, LinearAlgebra.dot, LinearAlgebra.normalize, LinearAlgebra.norm

export 
    Vector3,
    Matrix3,
    dot,
    cross,
    cross2,
    slerp,
    outer,
    inv,
    solve

include("Vector3.jl")
include("Matrix3.jl")


end 
    
