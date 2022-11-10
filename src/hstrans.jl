#=
hs变换
=#

using LinearAlgebra

"""
进行trotter分解
ret = exp(-Δτ*hk_Ltrot)exp(-Δτ*hu_Ltrot)...exp(-Δτ*hk_1)exp(-Δτ*hu_1)
顺序从右向左
"""
function trotter_ρ(hk, hu, β, Ltrot)
    ret = Vector{Matrix{ComplexF64}}(undef, 2*Ltrot)
    Δτ = β / Ltrot
    expdthk = exp(-Δτ*hk)
    expdthu = exp(-Δτ*hu)
    for tau in Base.OneTo(Ltrot)
        ret[2tau-1] = copy(expdthk)
        ret[2tau] = copy(expdthu)
    end
    return ret
end


"""
进行HS分解
创建一个exp(-ΔτU)
在Tr_{F} 这个费米子求trace的里面 0.5Tr_{F}(...) != Tr_{F}(0.5...)
所以常数权重单独拿出来
"""
@generated function hsdecomposite(Nsite, usite, Δτ)
    #γ::Vector{ComplexF64} = [1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    η::Vector{ComplexF64} = [-√(6+2√6), -√(6-2√6), √(6-2√6), √(6+2√6)]
    return quote
        ret = Array{Matrix{ComplexF64}, 2}(undef, Nsite, 4)
        for sidx in Base.OneTo(Nsite); for aidx in Base.OneTo(4)
            op = zeros(ComplexF64, Nsite*2, Nsite*2)
            op[sidx, sidx] = 1.0
            op[sidx+Nsite, sidx+Nsite] = 1.0
            #ret[sidx, aidx] = 0.25*$(γ)[aidx]*exp(-√(Δτ*usite[sidx])*$(η)[aidx]*op)
            ret[sidx, aidx] = exp(sqrt(Complex(-Δτ*usite[sidx]))*$(η)[aidx]*op)
        end; end
        ret
    end
end


"""
迭代所有hs configure
"""
struct HSIter
    #cfg :: Vector{Int, 2}
    ntau :: Int
    nsite :: Int
    length :: Int
    γ ::Vector{ComplexF64}
    HSIter(ntau, nsite) = new(
        #zeros(ntau, nsite),
        ntau,
        nsite,
        4^(ntau*nsite),
        [1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    )
end


Base.eltype(::HSIter) = Vector{Int, 2}
Base.length(hsit::HSIter) = hsit.length

function Base.iterate(hsit::HSIter, (count, cfg, wgt)=(0, missing, missing))
    if count >= hsit.length
        return nothing
    end 
    if count == 0
        cfg = ones(Int64, hsit.ntau, hsit.nsite)
        wgt = (hsit.γ[1]/4)^(hsit.ntau*hsit.nsite)
        return (cfg, wgt), (1, cfg, wgt)
    end
    for tidx in Base.OneTo(hsit.ntau); for sidx in Base.OneTo(hsit.nsite)
        if cfg[tidx, sidx] < 4
            cfg[tidx, sidx] += 1
            wgt = wgt * hsit.γ[cfg[tidx, sidx]] / hsit.γ[cfg[tidx, sidx]-1]
            @goto endloop
        else
            cfg[tidx, sidx] = 1
            wgt = wgt * hsit.γ[1] / hsit.γ[4]
        end
    end; end
    @label endloop
    return (cfg, wgt), (count+1, cfg, wgt)
end


"""
进行hs的逆变换，从
0.25*∑_{i=1,2,3,4} w_{i} e^{sqrt(ΔτU)ηA}
"""
function hscomposite(wgts, mats)
    throw(error("can not implement"))
end


"""
随机产生一个HS
"""
function random_configuration(Ltrot, Nsite)
    hsc = ceil.(Int64, 4rand(Ltrot, Nsite))
    return hsc
end


"""
B_l = e^{-ΔτK} e^{-ΔτV(l)}
"""
@generated function Bmat_τ(edk, vops, cffs, vhss)
    η::Vector{ComplexF64} = [-√(6+2√6), -√(6-2√6), √(6-2√6), √(6+2√6)]
    return quote
        bmat = edk
        vmat = zeros(size(vops[1]))
        for (op, cf, hs) in zip(vops, cffs, vhss)
            #println(hs)
            vmat = vmat + cf*$(η)[hs]*op
        end
        #println(vmat)
        bmat = bmat * exp(vmat)
    end
end


"""
从cfg计算Bseq
"""
function Bmarr(cfg::Matrix{T}, hub::GeneralHubbard) where T <: Integer
    Bseq = []
    siz = size(cfg)
    @assert length(hub.Hv) == length(hub.Hg) == siz[2]
    for lidx in Base.OneTo(siz[1])
        push!(Bseq, Bmat_τ(hub.e_dτHk, hub.Hv, hub.Hg, cfg[lidx, :]))
    end
    return Bseq
end
