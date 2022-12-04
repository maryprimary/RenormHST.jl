#=
hubbard相关的内容
=#


"""
构建一个hubbard模型的哈密顿
"""
hubbard_hamiltonian(hs, bonds, usite) = hubbard_hk(hs, bonds) + hubbard_hu(hs, usite)


"""
hubbard模型中的自由部分
"""
function hubbard_hk(hs, bonds)
    ham = zeros(ComplexF64, length(hs), length(hs))
    for idx1 in Base.OneTo(length(hs)); for idx2 in Base.OneTo(idx1)
        stl = hs[idx1]
        str = hs[idx2]
        for bnd in bonds
            if stl[bnd[1]] == '0' || str[bnd[2]] == '0'
                continue
            end
            lst = stl[1:bnd[1]-1]*"0"
            lsgn = 1
            for idx in 1:1:(bnd[1]-1)
                if stl[idx] == '1'
                    lsgn = lsgn * (-1)
                end
            end
            if bnd[1] < length(stl)
                lst = lst*stl[bnd[1]+1:end]
            end
            rst = str[1:bnd[2]-1]*"0"
            rsgn = 1
            for idx in 1:1:(bnd[2]-1)
                if str[idx] == '1'
                    rsgn = rsgn * (-1)
                end
            end
            if bnd[2] < length(str)
                rst = rst*str[bnd[2]+1:end]
            end
            if lst == rst
                ham[idx1, idx2] -= 1.0 * rsgn * lsgn
                ham[idx2, idx1] -= 1.0 * rsgn * lsgn
            end
        end
    end; end
    return ham
end


"""
hubbard模型的相互作用部分
"""
function hubbard_hu(hs, usite)
    ham = zeros(ComplexF64, length(hs), length(hs))
    for idx1 in Base.OneTo(length(hs))
        sta = hs[idx1]
        nsite = length(usite)
        for upair in enumerate(usite)
            if sta[upair[1]] == '1' && sta[upair[1]+nsite] == '1'
                ham[idx1, idx1] += upair[2] 
            end
        end
    end
    return ham
end


"""
通用的记录模型的方法
"""
struct GeneralHubbard
    Hk :: Matrix{ComplexF64}
    Hg :: Vector{ComplexF64}
    Hv :: Vector{Matrix{ComplexF64}}
    Δτ :: Float64
    e_dτHk :: Matrix{ComplexF64}
end


"""
产生一个
"""
function hubbard_onsite(hk, usite, Δτ)
    @assert size(hk)[1] == 2*length(usite)
    nsite = length(usite)
    cffs = [sqrt(Complex(-Δτ*u)) for u in usite]
    ops::Vector{Matrix{ComplexF64}} = []
    for idx in Base.OneTo(nsite)
        op = zeros(2*nsite, 2*nsite)
        op[idx, idx] = 1.0
        op[idx+nsite, idx+nsite] = 1.0
        push!(ops, op)
    end
    return GeneralHubbard(
        hk, cffs, ops, Δτ, exp(-Δτ*hk)
    )
end

