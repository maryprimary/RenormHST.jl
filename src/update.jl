#=
有关更新的算法
=#

using Random
using LinearAlgebra


mutable struct __updt_ctrl
    stblz :: Int64
end

update_control = __updt_ctrl(10)



"""
计算B的乘积
"""
function Bprod(Bseq)
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
翻转一个格子
"""
macro flip(cfg, cidx)
    return quote
        orihs = $(esc(cfg))[$(esc(cidx))]
        prop = ceil(Int64, 3rand())
        mod(orihs+prop-1, 4) + 1
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
