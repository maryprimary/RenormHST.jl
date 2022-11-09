using Test

include("../src/RenormHST.jl")
using ..RenormHST

using LinearAlgebra



@testset "Bmat" begin
    Bseq = []
    hk = [0 1; 1 0]
    Δτ = 0.1
    edk = exp(-Δτ*hk)
    hscfg = random_configuration(2, 1)
    println(hscfg)
    #
    v1 = [1 0; 0 1]
    vop = [v1]
    cfs = [sqrt(Complex(-Δτ*1.6))]
    #
    for lidx in Base.OneTo(2)
        push!(Bseq, Bmat_τ(edk, vop, cfs, hscfg[lidx, :]))
    end
    println(Bseq)
    println(RenormHST.update_control)
    RenormHST.update_control.stblz = 5
    println(RenormHST.update_control)
    println(@sgn 1 Bseq)
end


@testset "Z" begin
    #3.230381677739212 - 2.7755575615628914e-17im
    Δτ = 0.1
    v1 = [1 0; 0 1]
    vop = [v1]
    cfs = [sqrt(Complex(-Δτ*1.6))]
    hk = [0 0; 0 0]
    edk = exp(-Δτ*hk)
    #
    summ = 0.
    for (cfg, wgt) in HSIter(1, 1)
        Bseq = []
        for lidx in Base.OneTo(1)
            push!(Bseq, Bmat_τ(edk, vop, cfs, cfg[lidx, :]))
        end
        println(Bseq[1])
        summ = summ + wgt*(@sgn 1 Bseq)
    end
    println(summ)
end


@testset "Z2" begin
    #3.23145366768717 + 2.47198095326695e-17im
    #第二个
    Δτ = 0.05
    v1 = [1 0; 0 1]
    vop = [v1]
    cfs = [sqrt(Complex(-Δτ*1.6))]
    hk = [0 0; 0 0]
    edk = exp(-Δτ*hk)
    #
    summ = 0.
    for (cfg, wgt) in HSIter(2, 1)
        Bseq = []
        println(cfg, " ", cfg[:, 1])
        for lidx in Base.OneTo(2)
            push!(Bseq, Bmat_τ(edk, vop, cfs, cfg[lidx, :]))
        end
        println(Bseq[1], " ", wgt, " ", (@sgn 1 Bseq), " ", cfg)
        summ = summ + wgt*(@sgn 1 Bseq)
    end
    println(summ)
end


@testset "Bprod" begin
    Bseqt = []
    for idx in Base.OneTo(15)
        push!(Bseqt, rand(8, 8))
    end
    bpr = RenormHST.Bprod(Bseqt)
    bprudv = RenormHST.BprodUDV(Bseqt)
    @test all(isapprox.(bpr, bprudv, rtol=0.0001))
    println(RenormHST.update_control.stblz)
    RenormHST.update_control.stblz = 11
    bpr = RenormHST.Bprod(Bseqt)
    bprudv = RenormHST.BprodUDV(Bseqt)
    @test all(isapprox.(bpr, bprudv, rtol=0.0001))
end

