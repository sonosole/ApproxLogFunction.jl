using Test
using ApproxLogFunction


@testset "checking error requirement" begin
    N = 32;
    for type in [Float16, Float32, Float64],
        base in [1.997, π, 10],
        Δout in [0.2, 0.1, 0.05]
        alog = Approxlog(base, abserror=Δout, dtype=type);
        k = type(17)
        x = rand(type, N) .* k;
        y = alog.(x)
        z = log.(base, x)
        @test sum(y - z)/N ≤ Δout
    end
end


