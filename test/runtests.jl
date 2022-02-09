using Test
using LinearAlgebra
import Pfapack77: ltl_decomp!, ltl_inv

zrtest(val, atol, label) = begin
    iszr = ≈(val, 0.0, atol=atol)
    if !iszr
        @info "`$label` test failed. Amplitude $val"
    end
    return iszr
end

@testset "Pfaffian & Inv: Simple" begin

A = LowerTriangular(rand(400, 400)) ./ 3.5
B = A - A'

pf_amp = det(B'B)^0.25

L, T, iPrm = ltl_decomp!(B);
pf = prod(T[1:2:end])
@info pf
@info pf_amp
@test zrtest((abs(pf) - pf_amp)/pf_amp, 1e-10, "11_Pfaffian")

Bi = ltl_inv(L, T, iPrm)
@test zrtest(reduce(max, abs.(Bi * (A - A') - I)), 1e-8, "12_Inverse")

end
