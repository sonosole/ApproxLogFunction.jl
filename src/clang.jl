function csrc(funcname::String, f::Approxlog{I,F}) where {I,F}
    four  = "    "
    fmax  = floatmax(F)
    ftype = nothing
    itype = nothing
    logᵦ = f.logtable
    N    = length(logᵦ)
    code = "#include <stdint.h>\n\n"
    if I <: Int32
        ftype = "float"
        itype = "int32_t"
        code *= "static float table[$N] = \n{\n"
    end
    if I <: Int64
        ftype = "double"
        itype = "int64_t"
        code *= "static double table[$N] = \n{\n"
    end
    for i = 1 : N - 1
        code *= "$four$(logᵦ[i]),\n"
    end;code *= "$four$(logᵦ[N])\n};\n\n"

    fracbits  = f.fracbits
    expobias  = f.expobias
    fracmask  = f.fracmask
    shiftbits = f.shiftbits
    logbase2  = f.logbase2
    code *= """
    $ftype $funcname($ftype x)
    {
        if(x <= 0){
            return -$fmax;
        }
        $itype z = *( ($itype*) (&x) );
        $itype e = (z >> $fracbits) - $expobias;
        $itype i = (z & $fracmask) >> $shiftbits;
        return e * $logbase2 + table[i];
    }
    """
    return code
end


function hsrc(filename::String, funcname::String, f::Approxlog{I,F}) where {I,F}
    ftype = nothing
    if I <: Int32
        ftype = "float"
    end
    if I <: Int64
        ftype = "double"
    end

    FILENAME = uppercase(filename)
    head = """
    #ifndef $(FILENAME)_H
    #define $(FILENAME)_H

    #include "$(filename).h"

    $ftype $funcname($ftype x);

    #endif
    """
    return head
end


function toclang(dstdir::String, filename::String, funcname::String, f::Approxlog{I,F}) where {I,F}
    if !isdir(dstdir)
        mkdir(dstdir)
    end

    nowdir = pwd()
    cd(dstdir)

    str = hsrc(filename, funcname, f)
    dst = touch(filename * ".h")
    open(dst, "w") do io
        write(io, str)
    end

    str = csrc(funcname, f)
    dst = touch(filename * ".c")
    open(dst, "w") do io
        write(io, str)
    end

    cd(nowdir)
    return nothing
end


