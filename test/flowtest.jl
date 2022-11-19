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


#=
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
    println(sgn/ltm, " ", val/ltm, " ", val/sgn)
end
=#


#=
@testset "testshell" begin
    #shells0 = obser_shells([0.0])
    #shells1, shells2 = omega_shells(shells0)
    shells0 = [SignShell(-Inf, missing, 0.0), SignShell(0.0, missing, Inf)]
    shells1 = [SignShell(-Inf, 0.0, Inf)]
    println(shells0)
    println(shells1)
    #println(shells2)
    #
    function gfunc(Bseq)
        siz = size(Bseq[1])
        grf = inv(Diagonal(ones(siz[1])) + RenormHST.BprodUDV(Bseq))
        return grf
    end
    #
    Random.seed!(21174)
    hk = zeros(6, 6)
    Δτ = 0.01
    ohb = hubbard_onsite(hk, [1.6, 1.6, 1.6], Δτ)
    inicfg = random_configuration(10, 3)
    ltm = 10000
    #
    nums, dens = omegaflow(ltm, ohb, inicfg, shells1)
    sgn, val = obserflow(ltm, ohb, inicfg, shells0, gfunc)
    grfsu = val[1][2] / sgn[2]
    grfsd = val[1][1] / sgn[1]
    println(grfsu)
    println(grfsd)
    omega = nums[1] / dens[1]
    obser = (omega*grfsd + grfsu) / (1+omega)
    println(obser)
end
=#


@testset "shell" begin
    shells0 = obser_shells([-1e40, -1e30, -1e13, 0.0])
    shells1, shells2 = omega_shells(shells0)
    println(shells0)
    println(shells1)
    println(shells2)
    #
    function rho(grf)
        den = sum([grf[idx, idx] for idx in Base.OneTo(6)])
        return den / 6
    end
    #
    Random.seed!(21174)
    hk = zeros(6, 6)
    hk[1, 2] = -1.0
    hk[2, 1] = -1.0
    hk[4, 5] = -1.0
    hk[5, 4] = -1.0
    Δτ = 0.01
    ohb = hubbard_onsite(hk, [1.6, 1.6, 1.6], Δτ)
    inicfg = random_configuration(10, 3)
    ltm = 10000
    #
    oshell = Vector{SignShell}(undef, length(shells1)+length(shells2))
    omegas = Vector{ComplexF64}(undef, length(shells1)+length(shells2))
    #
    nums, dens = omegaflow(ltm, ohb, inicfg, shells1)
    for sidx in Base.OneTo(length(shells1))
        oshell[2*(sidx-1)+1] = shells1[sidx]
        omegas[2*(sidx-1)+1] = nums[sidx] / dens[sidx]
    end
    println(nums, " ", dens)
    nums, dens = omegaflow(ltm, ohb, inicfg, shells2)
    for sidx in Base.OneTo(length(shells2))
        oshell[2*sidx] = shells2[sidx]
        omegas[2*sidx] = nums[sidx] / dens[sidx]
    end   
    println(nums, " ", dens)
    println(oshell, " ", omegas)
    mshell, momega = omega_steps(oshell, omegas)
    println(mshell, momega)
    #
    sgn, grf, val = obserflow(ltm, ohb, inicfg, shells0, rho)
    dens = Vector{ComplexF64}(undef, length(shells0))
    for sidx in Base.OneTo(length(shells0))
        if isassigned(sgn, sidx) && isassigned(val[1], sidx)
            dens[sidx] = val[1][sidx] / sgn[sidx]
        end
    end
    mshell2, mobser = obser_steps(shells0, dens, momega)
    println(mshell2)
    for idx in 1:1:length(shells0)
        println(mobser[end-idx+1])
    end
end

#=
@testset "shell diff" begin
    #
    function rho(grf)
        den = sum([grf[idx, idx] for idx in Base.OneTo(6)])
        return den / 6
    end
    #
    Random.seed!(21174)
    hk = zeros(6, 6)
    Δτ = 0.01
    ohb = hubbard_onsite(hk, [1.6, 1.6, 1.6], Δτ)
    inicfg = random_configuration(10, 3)
    ltm = 10000
    #
    shells0 = obser_shells([-1e13, 0.0])
    sgn, grf, val = obserflow(ltm, ohb, inicfg, shells0, rho)
    grfs = Vector{ComplexF64}(undef, length(shells0))
    for sidx in Base.OneTo(length(shells0))
        if isassigned(sgn, sidx) && isassigned(val[1], sidx)
            grfs[sidx] = val[1][sidx] / sgn[sidx]
        end
    end
    println(grfs[end])
    #
    Random.seed!(21174)
    hk = zeros(6, 6)
    Δτ = 0.01
    ohb = hubbard_onsite(hk, [1.6, 1.6, 1.6], Δτ)
    inicfg = random_configuration(10, 3)
    ltm = 10000
    #
    shells0 = obser_shells([-1e14, -1e13, 0.0])
    sgn, grf, val = obserflow(ltm, ohb, inicfg, shells0, rho)
    println(sgn, val)
    grfs = Vector{ComplexF64}(undef, length(shells0))
    for sidx in Base.OneTo(length(shells0))
        if isassigned(sgn, sidx) && isassigned(val[1], sidx)
            grfs[sidx] = val[1][sidx] / sgn[sidx]
        end
    end
    println(grfs[end])
end
=#


#=
@testset "shell 3 step" begin
    shells0 = obser_shells([-1e13, -1e12, 0.0])
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
    Random.seed!(21174)
    hk = zeros(6, 6)
    Δτ = 0.01
    ohb = hubbard_onsite(hk, [1.6, 1.6, 1.6], Δτ)
    inicfg = random_configuration(10, 3)
    ltm = 10000
    #
    oshell = Vector{SignShell}(undef, length(shells1)+length(shells2))
    omegas = Vector{ComplexF64}(undef, length(shells1)+length(shells2))
    #
    nums, dens = omegaflow(ltm, ohb, inicfg, shells1)
    for sidx in Base.OneTo(length(shells1))
        oshell[2*(sidx-1)+1] = shells1[sidx]
        omegas[2*(sidx-1)+1] = nums[sidx] / dens[sidx]
    end
    println(nums, " ", dens)
    nums, dens = omegaflow(ltm, ohb, inicfg, shells2)
    for sidx in Base.OneTo(length(shells2))
        oshell[2*sidx] = shells2[sidx]
        omegas[2*sidx] = nums[sidx] / dens[sidx]
    end   
    println(nums, " ", dens)
    println(oshell, " ", omegas)
    mshell, momega = omega_steps(oshell, omegas)
    println(mshell, momega)
    #
    sgn, val = obserflow(ltm, ohb, inicfg, shells0, gfunc)
    grfs = Vector{Vector{ComplexF64}}(undef, length(shells0))
    for sidx in Base.OneTo(length(shells0))
        if isassigned(sgn, sidx) && isassigned(val[1], sidx)
            grfs[sidx] = diag(val[1][sidx]) / sgn[sidx]
        end
    end
    mshell2, mobser = obser_steps(shells0, grfs, momega)
    println(mshell2)
    println(mobser[4])
    println(mobser[3])
    println(mobser[2])
    println(mobser[1])
end
=#

