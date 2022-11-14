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
    omega = nums[1] / dens[1]
    obser = (omega*grfsd + grfsu) / (1+omega)
    println(obser)
end
=#



@testset "shell" begin
    shells0 = obser_shells([-1e5, 0])
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
    println(oshell)
    #
    nums, dens = omegaflow(1, ohb, inicfg, shells1)
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
    #
    sgn, val = obserflow(ltm, ohb, inicfg, shells0, gfunc)
    grfs = Vector{Vector{ComplexF64}}(undef, length(shells0))
    for sidx in Base.OneTo(length(shells0))
        if isassigned(sgn, sidx) && isassigned(val[1], sidx)
            grfs[sidx] = diag(val[1][sidx]) / sgn[sidx]
        end
    end
    #println(shells0[1], " ", sgn[1], " ", grfs[1])
    println(shells0[2], " ", sgn[2], " ", grfs[2])
    println(shells0[3], " ", sgn[3], " ", grfs[3])
    #
    obser = grfs[3]
    obser = (omegas[2]*grfs[2] + grfs[3]) / (1 + omegas[2])
    println(SignShell(oshell[2].lower, oshell[2].middl, Inf), obser, omegas[2])
    return
    omega = omegas[2]*omegas[1] / (1+omegas[1])
    obser = (omega*grfs[1] + obser) / (1 + omega)
    println(SignShell(oshell[1].lower, oshell[1].middl, Inf), obser, omega)
    return
    #println(sgn, val)
    #for (sidx, sh) in enumerate(shells0)
    #    println(sh, " ", sgn[sidx], " ", val[1][sidx])
    #    println(grfs[sidx])
    #end
    #
    nowshell = shells0[end]
    nowobser = grfs[end]
    nowomega = omegas[end]
    stepshell = oshell[end]
    fordshell = oshell[end-1]
    println(nowshell, "->", nowobser)
    for lidx in Base.OneTo(length(oshell)-2)
        println(nowshell, " ", stepshell, " ", fordshell)
        ob3, om3 = omega_step(nowobser, nowomega, grfs[end-lidx], omegas[end-lidx])
        newshell = SignShell(shells0[end-lidx].lower, shells0[end-lidx].upper, nowshell.upper)
        println(newshell, "->", ob3, " ", om3) 
        nowobser = ob3
        nowomega = om3
        nowshell = newshell
        stepshell = SignShell(oshell[end-lidx].lower, nowshell.lower, nowshell.upper)
        fordshell = oshell[end-lidx-1] 
    end
    println(nowshell, " ", stepshell, " ", fordshell)
    ob3, om3 = omega_step(nowobser, nowomega, grfs[end-length(oshell)+1], omegas[end-length(oshell)+1])
    newshell = SignShell(oshell[2].lower, oshell[2].upper, nowshell.upper)
    println(newshell, "->", ob3, " ", om3)
    nowobser = ob3
    nowomega = om3
    nowshell = newshell
    stepshell = SignShell(fordshell.lower, fordshell.middl, nowshell.upper)
    nowobser = (grfs[1]*nowomega + nowobser) / (1 + nowomega)
    println(nowomega, grfs[1], length(shells0)-length(oshell))
    println(stepshell, "->", nowobser)
end


