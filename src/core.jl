@inline function f2i(::Type{F}) where F <: AbstractFloat
    (F <: Float16) && return Int16
    (F <: Float32) && return Int32
    (F <: Float64) && return Int64
end


struct Approxlog{XInt, XFloat}
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


"""
    Approxlog(base::Real; dtype::Type{<:AbstractFloat}, abserror::Real)
A lookup table based method to calculate `y` = log(base, `x`) with controllable error

+ `base` is the radix of logarithm operation
+ `dtype` is the data type of input `x` and output `y`
+ `abserror` is the required absolute error of the output `y`

# Example
```julia
alog₂  = Approxlog(2, abserror=0.1, dtype=Float32);
input  = 5.20f2;
output = alog₂(input);
```
"""
function Approxlog(base::Real; dtype::Type{F}, abserror::Real) where F <: AbstractFloat
    @assert base > 0 "base > 0 shall be met, but got $base"
    @assert base ≠ 1 "base ≠ 1 shall be met, but got $base"

    if abs(1-base) < 1e-2
        @warn("the base of log is too much closing to 1")
    end

    I = f2i(F)  # the right integer type corresponding to float type F
    l = one(I)  # alias for integer 1

    inputerror   = abserror * abs(log(base))
    usedfracbits = floor(I, log2(nextpow(2, 1 / inputerror)))
    tablelength  = l << usedfracbits
    
    if usedfracbits > 14
        kbs = 2^usedfracbits * sizeof(F) / 1024
        @warn("the table buffer would occupy $kbs KiB")
    end

    Δx = F(1  / tablelength)     # the actual input error
    Δy = F(Δx / abs(log(base)))  # the actual output error

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
    expobias = ( (l << (expobits-1)) - l )  # exponent part's bias number for float type F
    fracmask = ( (l << (fracbits  )) - l )  # to set non-fracbits 0
    
    shiftbits = fracbits - usedfracbits
    @assert fracbits > usedfracbits "using too much bits for fraction, we only have $fracbits for $F, but you want $usedfracbits"

    table = Vector{F}(undef, tablelength)
    for i = 1 : tablelength
        table[i] = log(base, 1 + (i-1) * Δx)
    end

    return Approxlog{I,F}(expobits,
                          fracbits,
                          expobias,
                          fracmask,
                          shiftbits,
                          base, F(log(base,2)), Δy, table)
end


function (f::Approxlog{I,F})(x::F) where {I,F}
    x <  0 && error("input shall be greater than 0")
    x == 0 && return -floatmax(F)
    xint = reinterpret(I, x)
    e = (xint >> f.fracbits) - f.expobias
    i = (xint & f.fracmask) >> f.shiftbits
    return @inbounds e * f.logbase2 + f.logtable[i+1]
end


