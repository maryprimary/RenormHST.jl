#=
有关更新的算法
=#

using Random
using LinearAlgebra


mutable struct __updt_ctrl
    stblz :: Int64
    Bpdre :: Int64
end

update_control = __updt_ctrl(10, 1)



"""
计算B的乘积
"""
function BprodS(Bseq)
    ret = Bseq[1]
    for lidx in Base.OneTo(length(Bseq)-1)
        ret = Bseq[lidx+1] * ret
    end
    return ret
end



"""
包含stablize的B乘积
"""
function BprodUDV(Bseq)
    siz = size(Bseq[1])
    lst = 0
    #左侧积累U，右侧积累R
    U = Diagonal(ones(siz[1]))
    D = Diagonal(ones(siz[1]))
    R = Diagonal(ones(siz[1]))
    while lst < length(Bseq)
        remain = min(update_control.stblz, length(Bseq)-lst)
        for lidx in Base.OneTo(remain)
            U = Bseq[lst+lidx] * U
        end
        lst = lst+remain
        f = svd(U*D)
        U = f.U
        D = Diagonal(f.S)
        R = f.Vt * R
    end
    return U*D*R
end


"""
获取符号
A = wgt * ∏_{l} B_{l}
"""
macro sgn(wgt, Bseq)
    return quote
        _bseq = $(esc(Bseq))
        siz = size(_bseq[1])
        ind = Diagonal(ones(siz[1]))
        $(esc(wgt))*det(ind + BprodUDV(_bseq))
    end
end

"""
获取符号
"""
macro sgn2(wgt, Bprod)
    quote
        siz = size($(esc(Bprod)))
        ind = Diagonal(ones(siz[1]))
        $(esc(wgt))*det(ind + $(esc(Bprod)))
    end
end



"""
翻转一个格子
"""
macro flip(cfg, cidx)
    return quote
        orihs = $(esc(cfg))[$(esc(cidx))]
        prop = ceil(Int64, 3rand())
        mod(orihs+prop-1, 4) + 1
    end
end


"""
将Bseq推进一个，Bseq=[B(1),B(2),...], 变成Bseq=[B(2)...B(1)]。
需注意A=B(Ltrot)...B(2)B(1)数组和乘积是反过来的
"""
macro uptau_Bseq!(Bseq)
    quote
        #Bseq2 = Vector{Matrix{ComplexF64}}(undef, length($(esc(Bseq))))
        fl = $(esc(Bseq))[1]
        $(esc(Bseq))[1:end-1] = $(esc(Bseq))[2:end]
        $(esc(Bseq))[end] = fl
        $(esc(Bseq))
    end
end


"""
将Bprod推进一个
A=B(Ltrot)...B(2)B(1) -> A=B(1)B(Ltrot)...B(2)
"""
macro uptau_Bprod!(Bprod, Bseq)
    quote
        fl = $(esc(Bseq))[1]
        invfl = inv(fl)
        _bp = $(esc(Bprod))
        _bp = fl * _bp * invfl
        _bp
    end
end


"""
重新计算Bseq的第一个
"""
macro hspr_Bseq(Bseq, ohb, hsf)
    quote
        hub = $(esc(ohb))
        Bseq2 = Vector{Matrix{ComplexF64}}(undef, length($(esc(Bseq))))
        Bseq2[1] = Bmat_τ(hub.e_dτHk, hub.Hv, hub.Hg, $(esc(hsf)))
        Bseq2[2:end] = $(esc(Bseq))[2:end]
        Bseq2
    end
end


"""
更新一个Bmat
"""
macro hspr_Bmat(Bseq, ohb, site, hsnew, hsold)
    η::Vector{ComplexF64} = [-√(6+2√6), -√(6-2√6), √(6-2√6), √(6+2√6)]
    quote
        hub = $(esc(ohb))
        bmat = $(esc(Bseq))[1]
        vmat = hub.Hg[$(esc(site))]*($(η)[$(esc(hsnew))]-$(η)[$(esc(hsold))])*hub.Hv[$(esc(site))]
        #NOTE:
        #在我们的标记下，e^(Hv)在Bmat的右侧，这里假设了所有的Hv是对易的，于是直接乘
        bmat = bmat * exp(vmat)
        Bseq2 = Vector{Matrix{ComplexF64}}(undef, length($(esc(Bseq))))
        Bseq2[2:end] = $(esc(Bseq))[2:end]
        Bseq2[1] = bmat
        Bseq2
    end
end


"""
更新整个Bprod, A = B(Ltrot)...B(1) 
"""
macro hspr_Bprod(Bprod, ohb, site, hsnew, hsold)
    η::Vector{ComplexF64} = [-√(6+2√6), -√(6-2√6), √(6-2√6), √(6+2√6)]
    quote
        hub = $(esc(ohb))
        vmat = hub.Hg[$(esc(site))]*($(η)[$(esc(hsnew))]-$(η)[$(esc(hsold))])*hub.Hv[$(esc(site))]
        #NOTE:
        #在我们的标记下，e^(Hv)在Bmat的右侧，这里假设了所有的Hv是对易的，于是直接乘
        _bp = $(esc(Bprod)) * exp(vmat)
        _bp
    end
end

#=
"""
迭代所有的格点
"""
struct HSSweep 
    rwgt :: Float64
    initcfg :: Matrix{Int}
    ntau :: Int
    nsite :: Int
    length :: Int
    γ ::Vector{ComplexF64}
    HSSweep(initc, seed) = new(
        1.0,
        initc,
        size(initc)[1],
        size(initc)[2],
        prod(size(initc)),
        [1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    )
end


Base.eltype(::HSSweep) = Vector{Int, 2}
Base.length(hsit::HSSweep) = hsit.length


function Base.iterate(hsit::HSSweep, (count, cfg, wgt)=(0, missing, missing))
    if count >= hsit.length
        return nothing
    end 
    if count == 0
        cfg = copy(hsit.initcfg)
        wgt = hsit.rwgt
        return (cfg, wgt), (1, cfg, wgt)
    end
    orihs = cfg[count]
    prop = ceil(Int64, 3rand())
    prop = mod(orihs+prop-1, 4) + 1
    
    return (cfg, wgt), (count+1, cfg, wgt)
end
=#

