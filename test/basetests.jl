using Test

include("../src/RenormHST.jl")
using ..RenormHST

using LinearAlgebra


@testset "Hilbert" begin
    println(fermi_hilbert(2))
    println(fermi_hilbert(3))
end


@testset "Hubbard" begin
    hs = fermi_hilbert(4)
    hh = hubbard_hamiltonian(hs, [(1, 2), (2, 1), (3, 4), (4, 3)], [0., 0.])
    #println(hh)
    #for idx1 in Base.OneTo(length(hs))
    #    for idx2 in Base.OneTo(idx1)
    #        stl = hs[idx1]
    #        str = hs[idx2]
    #        println(stl, " ", str, " ", hh[idx1, idx2], " ", hh[idx2, idx1])
    #    end
    #    sta = hs[idx1]
    #    println(sta, " ", hh[idx1, idx1])
    #end
    engr = eigvals(hh)
    println(engr)
end


@testset "trace" begin
    hs = fermi_hilbert(4)
    hh = hubbard_hamiltonian(hs, [(1, 2), (2, 1), (3, 4), (4, 3)], [0., 0.])
    β = 1.0
    ρ = exp(-β*hh)
    trace = tr(ρ)
    #
    ind = Diagonal(ones(4))
    hsp = tightbinding_hamiltonian(2, [(1, 2), (2, 1)], [0., 0.])
    tr2 = det(ind+exp(-β*hsp))
    println(trace, " ", tr2)
end


@testset "trotter" begin
    hs = fermi_hilbert(4)
    hk = hubbard_hk(hs, [(1, 2), (2, 1), (3, 4), (4, 3)])
    hu = hubbard_hu(hs, [2., 2.])
    hh = hk + hu
    #
    β = 1.0
    ρ = exp(-β*hh)
    trace = tr(ρ)
    #
    slices = trotter_ρ(hk, hu, β, 8)
    ind = Diagonal(ones(2))
    tr2 = tr(prod(slices))
    #
    slices2 = trotter_ρ(hk, hu, β, 80)
    ind = Diagonal(ones(2))
    tr3 = tr(prod(slices2))
    println(trace, " ", tr2, " ", tr3)
end


@testset "hst" begin
    Δτ = 0.1
    #
    hs = fermi_hilbert(2)
    hk = hubbard_hk(hs, [(1, 1), (2, 2)])
    hu = hubbard_hu(hs, [3.2])
    hh = -0.8*hk + hu
    ρ = exp(-Δτ*hh)
    trace = tr(ρ)
    #
    hst = hsdecomposite(1, [1.6], Δτ)
    println(trace)
    println(hst)
    println(size(hst))
    summ = 0.
    ind = Diagonal(ones(2))
    coef = [1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    for (cfg, wgt) in HSIter(1, 1)
        println(cfg, size(cfg))
        mat = Diagonal(ones(2))
        for tidx in Base.OneTo(1); for sidx in Base.OneTo(1)
            mat = mat*hst[sidx, cfg[1, sidx]]
        end; end
        println(mat)
        summ += wgt*det(ind + mat)
    end
    println(summ)
end

