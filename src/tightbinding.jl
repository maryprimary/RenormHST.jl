#=
紧束缚的哈密顿量
=#


"""
紧束缚哈密顿量
"""
tightbinding_hamiltonian(Nsite, bonds, usite) =
tightbinding_hk(Nsite, bonds) + tightbinding_hu(Nsite, usite)



"""
hopping
"""
function tightbinding_hk(Nsite, bonds)
    ham = zeros(ComplexF64, Nsite*2, Nsite*2)
    for bnd in bonds
        ham[bnd[1], bnd[2]] -= 1.0
        ham[bnd[1]+Nsite, bnd[2]+Nsite] -= 1.0
    end
    return ham
end


"""
interaction
"""
function tightbinding_hu(Nsite, usite)
    ham = zeros(ComplexF64, Nsite*2, Nsite*2)
    for upair in enumerate(usite)
        ham[upair[1], upair[1]] += upair[2]
        ham[upair[1]+Nsite, upair[1]] += upair[2]
    end
    return ham
end



