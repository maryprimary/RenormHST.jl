#=
测试flow
=#


using Test


include("../src/RenormHST.jl")
using ..RenormHST

using Random
using LinearAlgebra


@testset "exact 4 site" begin
    hs = fermi_hilbert(8)
    hu = 1.6
    println(length(hs))
    #-2
    hmu = hubbard_hk(hs, [(idx, idx) for idx in 1:1:8])
    #
    hk = hubbard_hk(hs, [(1, 2), (2, 1), (1, 3), (3, 1), (1, 4), (4, 1),
    (2, 3), (3, 2), (2, 4), (4, 2), (3, 4), (4, 3),
    (5, 6), (6, 5), (5, 7), (7, 5), (5, 8), (8, 5),
    (6, 7), (7, 6), (6, 8), (8, 6), (7, 8), (8, 7)])
    #hk1 = hubbard_hk(hs, [(1, 1)])
    #hk2 = hubbard_hk(hs, [(2, 2)])
    hu = hubbard_hu(hs, [8.0, 8.0, 8.0, 8.0])
    hh = hu + hk + 2.0*hmu
    rho = exp(-10.0*hh)
    den = zeros(8)
    for (cidx, cfg) in enumerate(hs)
        for idx in Base.OneTo(8)
            den[idx] += cfg[idx] == '1' ? rho[cidx, cidx] : 0
        end
    end
    println(tr(rho), den/tr(rho))
end


@testset "exact 6 site" begin
    hs = fermi_hilbert(12)
    #-2
    hmu = hubbard_hk(hs, [(idx, idx) for idx in 1:1:12])
    #
    hk = hubbard_hk(hs, [(1, 2), (1, 3), (1, 5), (1, 6), (2, 3), (2, 4),
    (2, 6), (2, 1), (3, 4), (3, 5), (3, 1), (3, 2), (4, 5), (4, 6),
    (4, 2), (4, 3), (5, 6), (5, 1), (5, 3), (5, 4),
    (6, 1), (6, 2), (6, 4), (6, 5),
    (7, 8), (7, 9), (7, 11), (7, 12), (8, 9), (8, 10),
    (8, 12), (8, 7), (9, 10), (9, 11), (9, 7), (9, 8), (10, 11), (10, 12),
    (10, 8), (10, 9), (11, 12), (11, 7), (11, 9), (11, 10),
    (12, 7), (12, 8), (12, 10), (12, 11)])
    #hk1 = hubbard_hk(hs, [(1, 1)])
    #hk2 = hubbard_hk(hs, [(2, 2)])
    hu = hubbard_hu(hs, [6.0, 6.0, 6.0, 6.0, 6.0, 6.0])
    hh = hu + hk + 1.5*hmu
    rho = exp(-6.0*hh)
    den = zeros(12)
    for (cidx, cfg) in enumerate(hs)
        for idx in Base.OneTo(12)
            den[idx] += cfg[idx] == '1' ? rho[cidx, cidx] : 0
        end
    end
    println(tr(rho), den/tr(rho))
end

#=
@testset "mc ord" begin
    #
    function rho(grf)
        den = sum([grf[idx, idx] for idx in Base.OneTo(6)])
        return den / 6
    end
    #
    hk = zeros(6, 6)
    Δτ = 0.1
    ohb = hubbard_onsite(hk, [1.6, 1.6, 1.6], Δτ)
    inicfg = random_configuration(5, 3)
    ltm = 100000
    sgn, grf, val = mcsum(ltm, ohb, inicfg, rho)
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
    function rho(grf)
        den = sum([grf[idx, idx] for idx in Base.OneTo(6)])
        return den / 6
    end
    #
    Random.seed!(21174)
    hk = zeros(6, 6)
    Δτ = 0.1
    ohb = hubbard_onsite(hk, [1.6, 1.6, 1.6], Δτ)
    inicfg = random_configuration(5, 3)
    ltm = 5000
    #
    nums, dens = omegaflow(ltm, ohb, inicfg, shells1)
    sgn, grf, val = obserflow(ltm, ohb, inicfg, shells0, rho)
    println(sgn, val)
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
    shells0 = obser_shells([-1e30])
    shells1, shells2 = omega_shells(shells0)
    println(shells0)
    println(shells1)
    println(shells2)
    #
    function rho(grf)
        den = sum([grf[idx, idx] for idx in Base.OneTo(8)])
        return den / 8
    end
    #
    Random.seed!(21251)
    hk = zeros(8, 8)
    #化学势
    #for idx in 1:1:8
    #    hk[idx, idx] = -3.2
    #end
    hbonds = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4),
    (3, 4), (5, 6), (5, 7), (5, 8), (6, 7), (6, 8), (7, 8)]
    #hopping
    for hop in hbonds
        hk[hop[1], hop[2]] += -1.0
        hk[hop[2], hop[1]] += -1.0
    end
    println(hk)
    println(eigvals(hk))
    Δτ = 0.1
    ohb = hubbard_onsite(hk, [-0.5, -0.5, -0.5, -0.5], Δτ)
    inicfg = random_configuration(40, 4)
    inicfg = warmup_configuration(10000, ohb, inicfg)
    ltm = 20000
    btm = 4
    #println(sgn_hist(ltm, ohb, inicfg, 500))
    #return
    #求所有的omega
    oshell = Vector{SignShell}(undef, length(shells1)+length(shells2))
    #omegas = Vector{ComplexF64}(undef, length(shells1)+length(shells2))
    #omegas_err = Vector{ComplexF64}(undef, length(shells1)+length(shells2))
    for _sidx in Base.OneTo(length(shells0)-1)
        sidx = length(shells0) - _sidx
        oshell[sidx] = SignShell(shells0[sidx].lower, shells0[sidx].upper, shells0[sidx+1].upper)
    end
    #
    sgn = []
    sgn_err = []
    sgn_bin = []
    omega = []
    omega_err = []
    omega_bin = []
    val = []
    val_err = []
    val_bin = []
    for bidx in Base.OneTo(btm)
        _sgn, _grf, _val = obserflow(ltm, ohb, inicfg, shells0, rho)
        println(_sgn)
        push!(sgn_bin, _sgn)
        push!(val_bin, _val[1] ./ _sgn)
        println(_val[1] ./ _sgn)
        _omg = Vector{ComplexF64}(undef, length(shells1)+length(shells2))
        for _sidx in Base.OneTo(length(shells0)-1)
            sidx = length(shells0) - _sidx
            _omg[sidx] = _sgn[sidx] / _sgn[sidx+1]
        end
        push!(omega_bin, _omg)
    end
    #整理平均值和误差
    for sidx in Base.OneTo(length(shells0))
        sgnb = [sgn_bin[bidx][sidx] for bidx in Base.OneTo(btm)]
        valb = [val_bin[bidx][sidx] for bidx in Base.OneTo(btm)]
        push!(sgn, sum(sgnb) / btm)
        push!(val, sum(valb) / btm)
        serr = (sgnb .- sgn[end]).^2
        push!(sgn_err, sqrt(sum(serr) / (btm-1)))
        verr = (valb .- val[end]).^2
        push!(val_err, sqrt(sum(verr) / (btm-1)))
    end
    for sidx in Base.OneTo(length(shells0)-1)
        omgb = [omega_bin[bidx][sidx] for bidx in Base.OneTo(btm)]
        push!(omega, sum(omgb) / btm)
        oerr = (omgb .- omega[end]).^2
        push!(omega_err, sqrt(sum(oerr) / (btm-1)))
    end
    println("shell中的符号 ", sgn, sgn_err)
    println("shell中的观测 ", val, val_err)
    println("shell中的流", omega, omega_err)
    println("约束Omega", oshell, " ", omega, omega_err)
    mshell, momega, momega_err = omega_steps(oshell, omega, omega_err)
    println("不约束的Omega", mshell, " ", momega, " ", momega_err)
    #sgn, grf, _val = obserflow(ltm, ohb, inicfg, shells0, rho)
    #dens = Vector{ComplexF64}(undef, length(shells0))
    #for sidx in Base.OneTo(length(shells0))
    #    #if isassigned(sgn, sidx) && isassigned(val[1], sidx)
    #        dens[sidx] = _val[1][sidx] / sgn[sidx]
    #    #end
    #end
    #println(dens)
    mshell2, mobser, mobser_err = obser_steps(shells0, val, momega, val_err, momega_err)
    println(mshell2)
    println("sgn ", sgn)
    for idx in 1:1:length(shells0)
        println(mobser[end-idx+1], " ", mobser_err[end-idx+1])
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

