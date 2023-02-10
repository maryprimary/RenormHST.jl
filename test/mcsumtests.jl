using Test


include("../src/RenormHST.jl")
using ..RenormHST

using Random
using LinearAlgebra


@testset "exact" begin
    hs = fermi_hilbert(2)
    hu = 1.6
    println(hs)
    den = [0.0, 0.0]
    Z = 0.
    for cfg in hs
        eng = 0.
        if cfg[1] == '1'
            eng += hu
        end
        if cfg[2] == '1'
            eng += hu
        end
        if cfg == "11"
            eng += 2*hu
        end
        ebh = exp(-0.1*eng)
        den[1] += cfg[1] == '1' ? ebh : 0
        den[2] += cfg[2] == '1' ? ebh : 0
        Z += ebh
    end
    println("den ", den/Z)
end


@testset "mcsum" begin
    function rho(grf)
        return 0.5*(grf[1, 1] + grf[2, 2])
    end
    Random.seed!(23441)
    hk = [0.0 0.0; 0.0 0.0]
    Δτ = 0.02
    ohb = hubbard_onsite(hk, [1.6], Δτ)
    println(ohb)
    inicfg = random_configuration(5, 1)
    inicfg = warmup_configuration(1000, ohb, inicfg)
    gamma = 0.25*[1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    wgtori = prod(map(x->gamma[x], inicfg))
    println(wgtori, " wgtori")
    ltm = 10000
    sgn, grf, val = mcsum(ltm, ohb, inicfg, rho)
    #println(16*val*wgtori/ltm)
    println(sgn/ltm, " ", grf/sgn, " ", val/sgn)
end

#=
@testset "mcflow" begin
    function gfunc(Bseq)
        siz = size(Bseq[1])
        grf = inv(Diagonal(ones(siz[1])) + RenormHST.BprodUDV(Bseq))
        return grf
    end
    hk = [0 0; 0 0]
    Δτ = 0.01
    ohb = hubbard_onsite(hk, [1.6], Δτ)
    println(ohb)
    inicfg = random_configuration(10, 1)
    gamma = 0.25*[1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    wgtori = prod(map(x->gamma[x], inicfg))
    ltm = 10000
    #shells
    shells::Vector{SignShell} = []
    push!(shells, SignShell(-Inf, missing, -1e8))
    push!(shells, SignShell(-1e8, missing, -1e7))
    push!(shells, SignShell(-1e7, missing, Inf))
    #
    sgn, val = obserflow(ltm, ohb, inicfg, shells, gfunc)
    #println(16*val*wgtori/ltm)
    println(sgn, " ", val)
    for (sidx, sh) in enumerate(shells)
        println(sh)
        println(sgn[sidx]/ltm)
        println(val[1])
        if isassigned(val[1], sidx)
            println(val[1][sidx]/ltm)
        end 
    end
end
=#

#=
@testset "mcflow2" begin
    function gfunc(Bseq)
        siz = size(Bseq[1])
        grf = inv(Diagonal(ones(siz[1])) + RenormHST.BprodUDV(Bseq))
        return grf
    end
    hk = [0 0; 0 0]
    Δτ = 0.01
    ohb = hubbard_onsite(hk, [1.6], Δτ)
    println(ohb)
    inicfg = random_configuration(10, 1)
    gamma = 0.25*[1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    wgtori = prod(map(x->gamma[x], inicfg))
    ltm = 10000
    #shells
    shells1::Vector{SignShell} = []
    push!(shells1, SignShell(-Inf, missing, -1e8))
    push!(shells1, SignShell(-1e8, missing, 0))
    push!(shells1, SignShell(0, missing, Inf))
    shells2::Vector{SignShell} = []
    for sidx in Base.OneTo(length(shells1)-1)
        push!(shells2, SignShell(shells1[sidx].lower, shells1[sidx].upper, shells1[sidx+1].upper))
    end
    println(shells1)
    println(shells2)
    #
    #sgn, val = obserflow(ltm, ohb, inicfg, shells1, gfunc)
    #println(sgn, val)
    #println(val[1][2] / sgn[2])
    #println(val[1][3] / sgn[3])
    Random.seed!(23124)
    nums, dens = omegaflow(ltm, ohb, inicfg, shells2)
    println(nums, dens)
    #nums, dens = obserflow(ltm, ohb, inicfg, shells, gfunc)
    ##println(16*val*wgtori/ltm)
    #println(sgn, " ", val)
    #for (sidx, sh) in enumerate(shells)
    #    println(sh)
    #    println(sgn[sidx]/ltm)
    #    println(val[1])
    #    if isassigned(val[1], sidx)
    #        println(val[1][sidx]/ltm)
    #    end 
    #end
end
=#
