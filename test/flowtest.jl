#=
测试flow
=#


using Test


include("../src/RenormHST.jl")
using ..RenormHST

using Random
using LinearAlgebra


@testset "exact 3 site" begin
    hs = fermi_hilbert(6)
    hu = 1.6
    println(length(hs))
    den = zeros(6)
    Z = 0.
    for cfg in hs
        eng = 0.
        for idx in Base.OneTo(6)
            if cfg[idx] == '1'
                eng += hu
            end
        end
        for idx in Base.OneTo(3)
            if cfg[idx] == '1' && cfg[idx+3] == '1'
                eng += 2*hu
            end
        end
        ebh = exp(-0.1*eng)
        for idx in Base.OneTo(6)
            den[idx] += cfg[idx] == '1' ? ebh : 0
        end
        Z += ebh
    end
    println("den ", den/Z)
end



@testset "mc ord" begin
    #
    function gfunc(Bseq)
        siz = size(Bseq[1])
        grf = inv(Diagonal(ones(siz[1])) + RenormHST.BprodUDV(Bseq))
        return grf
    end
    #
    hk = zeros(6, 6)
    Δτ = 0.01
    ohb = hubbard_onsite(hk, [1.6, 1.6, 1.6], Δτ)
    inicfg = random_configuration(10, 3)
    ltm = 10000
    sgn, val = mcsum(ltm, ohb, inicfg, gfunc)
    #println(16*val*wgtori/ltm)
    println(sgn/ltm, " ", val/ltm)
end


#=
@testset "shell" begin
    shells0 = obser_shells([-1e8, -1e7, -1e6, 0])
    shells1, shells2 = omega_shells(shells0)
    println(shells0)
    println(shells1)
    println(shells2)
    #
    function gfunc(Bseq)
        siz = size(Bseq[1])
        grf = inv(Diagonal(ones(siz[1])) + RenormHST.BprodUDV(Bseq))
        return grf
    end
    #
    hk = zeros(6, 6)
    Δτ = 0.01
    ohb = hubbard_onsite(hk, [1.6, 1.6, 1.6], Δτ)
    inicfg = random_configuration(10, 3)
    ltm = 10000
    #
    nums, dens = omegaflow(ltm, ohb, inicfg, shells1)
    println(nums, " ", dens)
    nums, dens = omegaflow(ltm, ohb, inicfg, shells2)
    println(nums, " ", dens)
    #
    sgn, val = obserflow(ltm, ohb, inicfg, shells0, gfunc)
    println(sgn, val)
end
=#

