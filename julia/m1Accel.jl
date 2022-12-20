using AppleAccelerateLinAlgWrapper
using AppleAccelerate
using DSP

#AppleAccelerate.@replaceBase round exp exp2 expm1 log log1p log2 log10 sin cos cospi tan asin acos atan cis sinh cosh tanh asinh acosh atanh sqrt copysign exponent abs rem ceil floor

import Base.*
*(A::Matrix{Float64}, B::Matrix{Float64}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{Float32}, B::Matrix{Float32}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{Float16}, B::Matrix{Float16}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{Bool}, B::Matrix{Bool}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{Int8}, B::Matrix{Int8}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{UInt8}, B::Matrix{UInt8}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{Int16}, B::Matrix{Int16}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{UInt16}, B::Matrix{UInt16}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{Int32}, B::Matrix{Int32}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{UInt32}, B::Matrix{UInt32}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{Int64}, B::Matrix{Int64}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{UInt64}, B::Matrix{UInt64}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{Int128}, B::Matrix{Int128}) = AppleAccelerateLinAlgWrapper.gemm(A, B)
*(A::Matrix{UInt128}, B::Matrix{UInt128}) = AppleAccelerateLinAlgWrapper.gemm(A, B)

#= import Base./
/(A::Matrix{Float64}, B::Matrix{Float64}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{Float32}, B::Matrix{Float32}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{Float16}, B::Matrix{Float16}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{Bool}, B::Matrix{Bool}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{Int8}, B::Matrix{Int8}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{UInt8}, B::Matrix{UInt8}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{Int16}, B::Matrix{Int16}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{UInt16}, B::Matrix{UInt16}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{Int32}, B::Matrix{Int32}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{UInt32}, B::Matrix{UInt32}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{Int64}, B::Matrix{Int64}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{UInt64}, B::Matrix{UInt64}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{Int128}, B::Matrix{Int128}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)
/(A::Matrix{UInt128}, B::Matrix{UInt128}) = AppleAccelerateLinAlgWrapper.ldiv(A, B)

import Base.\
\(A::Matrix{Float64}, B::Matrix{Float64}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{Float32}, B::Matrix{Float32}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{Float16}, B::Matrix{Float16}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{Bool}, B::Matrix{Bool}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{Int8}, B::Matrix{Int8}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{UInt8}, B::Matrix{UInt8}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{Int16}, B::Matrix{Int16}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{UInt16}, B::Matrix{UInt16}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{Int32}, B::Matrix{Int32}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{UInt32}, B::Matrix{UInt32}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{Int64}, B::Matrix{Int64}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{UInt64}, B::Matrix{UInt64}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{Int128}, B::Matrix{Int128}) = AppleAccelerateLinAlgWrapper.rdiv(A, B)
\(A::Matrix{UInt128}, B::Matrix{UInt128}) = AppleAccelerateLinAlgWrapper.rdiv(A, B) =#