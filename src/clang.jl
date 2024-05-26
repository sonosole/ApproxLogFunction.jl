function hcsrc(filename::String, funcname::String, f::Approxlog{I,F}) where {I,F}
    four  = "    "
    fmax  = floatmax(F)
    ftype = nothing
    itype = nothing
    logᵦ = f.logtable
    N    = length(logᵦ)

    dotc = """
    #include <stdint.h>
    #include "$(filename).h"\n
    """

    if I <: Int32
        ftype = "float"
        itype = "int32_t"
        dotc *= "static float table[$N] = \n{\n"
    end
    if I <: Int64
        ftype = "double"
        itype = "int64_t"
        dotc *= "static double table[$N] = \n{\n"
    end
    for i = 1 : N - 1
        dotc *= "$four$(logᵦ[i]),\n"
    end;dotc *= "$four$(logᵦ[N])\n};\n\n"

    fracbits  = f.fracbits
    expobias  = f.expobias
    fracmask  = f.fracmask
    shiftbits = f.shiftbits
    logbase2  = f.logbase2
    dotc *= """
    $ftype $funcname($ftype x)
    {
        if(x <= 0){
            return -$fmax;
        }
        $itype z = *( ($itype*) (&x) );
        $itype e = (z >> $fracbits) - $expobias;
        $itype i = (z & $fracmask) >> $shiftbits;
        return $logbase2 * e + table[i];
    }
    """

    FILENAME = uppercase(filename)
    doth = """
    #ifndef $(FILENAME)_H
    #define $(FILENAME)_H

    $ftype $funcname($ftype x);

    #endif
    """
    return doth, dotc
end


"""
    toclang(dstdir::String, filename::String, funcname::String, f::Approxlog)
generate C language `.h` and `.c` files for an `Approxlog` functor `f`.

+ `dstdir` is the directory to save the C language files (`dir` would be made if it doesn't exist before)
+ `filename` is the name of the generated `.c` and `.h` files
+ `funcname` is the name of the generated approximate function

# Example
```julia
alog₂ = Approxlog(2, abserror=0.12, dtype=Float32);
toclang("./cfolder/",  "approxlog", "alog", alog₂);
```
"""
function toclang(dstdir::String, filename::String, funcname::String, f::Approxlog{I,F}) where {I,F}
    if !isdir(dstdir)
        mkdir(dstdir)
    end

    nowdir = pwd()
    cd(dstdir)

    doth, dotc = hcsrc(filename, funcname, f)

    hdst = touch(filename * ".h")
    cdst = touch(filename * ".c")

    open(hdst, "w") do io
        write(io, doth)
    end

    open(cdst, "w") do io
        write(io, dotc)
    end

    cd(nowdir)
    return nothing
end


