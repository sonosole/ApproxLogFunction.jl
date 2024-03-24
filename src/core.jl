@inline function f2i(::Type{F}) where F <: AbstractFloat
    (F <: Float16) && return Int16
    (F <: Float32) && return Int32
    (F <: Float64) && return Int64
end


mutable struct Approxlog{XInt, XFloat}
    expobits  :: XInt
    fracbits  :: XInt
    expobias  :: XInt
    fracmask  :: XInt
    shiftbits :: XInt
    logbase   :: Real
    logbase2  :: XFloat
    maxerror  :: XFloat
    logtable  :: Vector{XFloat}
    function Approxlog{XInt,XFloat}(expobits  :: XInt,
                                    fracbits  :: XInt,
                                    expobias  :: XInt,
                                    fracmask  :: XInt,
                                    shiftbits :: XInt,
                                    logbase   :: Real,
                                    logbase2  :: XFloat,
                                    maxerror  :: XFloat,
                                    logtable  :: Vector{XFloat}) where {XInt <: Integer, XFloat <: AbstractFloat}

        new{XInt,XFloat}(expobits, fracbits, expobias, fracmask,
                 shiftbits, logbase, logbase2, maxerror, logtable)
    end
end


function Base.show(io::IO, ::MIME"text/plain", a::Approxlog{I,F}) where {I,F}
    logbase  = a.logbase
    maxerror = a.maxerror
    expobits = a.expobits
    fracbits = a.fracbits
    expobias = a.expobias
    fracmask = a.fracmask
    shiftbits = a.shiftbits
    println("Approxlog{$logbase, $F}:")
    println("  maxerror  = $maxerror")
    println("  signbits  = 1")
    println("  expobits  = $expobits")
    println("  fracbits  = $fracbits")
    println("  expobias  = $expobias")
    println("  fracmask  = $fracmask")
      print("  shiftbits = $shiftbits")
end


function table(f::Approxlog)
    return f.logtable
end


function maxerror(f::Approxlog)
    return f.maxerror
end


function Approxlog(base::Real; dtype::Type{F}, abserror::Real) where F <: AbstractFloat
    @assert base > 1 "base > 1 shall be met, but got $base"
    I  = f2i(F)
    l  = one(I)
    dx = abserror * log(base)

    npower = nextpow(2, 1 / dx)
    usedfracbits = floor(I, log2(npower))
    
    tsize = l << usedfracbits
    table = Vector{F}(undef, tsize)
    delta = 1 / tsize
    for i = 1 : tsize
        table[i] = log(base, 1 + (i-1) * delta)
    end

    maxerror = delta / log(base)
    expobits = nothing
    fracbits = nothing
    if F <: Float16
        expobits = I(5)
        fracbits = I(10)
    end
    if F <: Float32
        expobits = I(8)
        fracbits = I(23)
    end
    if F <: Float64
        expobits = I(11)
        fracbits = I(52)
    end

    expobias  = ( (l << (expobits-1)) - l )
    fracmask  = ( (l << (fracbits  )) - l )
    shiftbits = fracbits - usedfracbits
    return Approxlog{I,F}(expobits,
                          fracbits,
                          expobias,
                          fracmask,
                          shiftbits,
                          base,
                          F(log(base,2)),
                          F(maxerror),
                          table)
end


function (f::Approxlog{I,F})(x::F) where {I,F}
    x <  0 && error("input shall be greater than 0")
    x == 0 && return -floatmax(F)
    xint = reinterpret(I, x)
    e = (xint >> f.fracbits) - f.expobias
    i = (xint & f.fracmask) >> f.shiftbits
    return @inbounds e * f.logbase2 + f.logtable[i+1]
end


