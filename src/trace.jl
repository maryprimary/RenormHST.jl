#=
生成trace相关的功能
=#



"""
费米子Hilbert空间
"""
function fermi_hilbert(Nsite)
    function __gen_next(stas, ::Val{0})
        return stas
    end
    function __gen_next(stas, ::Val{idx}) where idx
        ret::Vector{String} = []
        for sta in stas
            push!(ret, sta)
            push!(ret, sta[1:idx-1]*"1"*sta[idx+1:end])
        end
        return __gen_next(ret, Val(idx-1))
    end
    return __gen_next([repeat("0", Nsite)], Val(Nsite))
end




