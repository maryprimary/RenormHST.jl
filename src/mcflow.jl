#=
进行mc的计算
=#


struct SignShell
    lower :: Float64
    middl :: Union{Float64, Missing}
    upper :: Float64
end


"""
推到合适的位置
"""
macro push_shell(ReSC, shells, arr, val)
    quote
        for (sidx, sh) in enumerate($(esc(shells)))
            if $(esc(ReSC)) < sh.upper && $(esc(ReSC)) >= sh.lower
                $(esc(arr))[sidx] = $(esc(val))
            end
        end
    end
end


"""
加到合适的位置
"""
macro add_shell(ReSC, shells, arr, val)
    quote
        for (sidx, sh) in enumerate($(esc(shells)))
            if $(esc(ReSC)) < sh.upper && $(esc(ReSC)) >= sh.lower
                $(esc(arr))[sidx] += $(esc(val))
            end
        end
    end
end


"""
找到合适的位置
"""
macro pull_shell(ReSC, shells, arr)
    quote
        val = missing
        for (sidx, sh) in enumerate($(esc(shells)))
            if $(esc(ReSC)) < sh.upper && $(esc(ReSC)) >= sh.lower
                val = $(esc(arr))[sidx]
            end
        end
        val
    end
end


"""
进行looptime次的mc step，可以作为一个bin
"""
function obserflow(looptime::T, 
    hub::GeneralHubbard, inicfg::Matrix{P}, shells::Vector{SignShell},
    valexprs...) where {T<:Integer, P<:Integer}
    #
    _γ = [1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    oldcfg = copy(inicfg)
    oldwgt = 1.0
    prcfg = copy(inicfg)
    Bseq = Bmarr(oldcfg, hub)
    oldSC = @sgn 1.0 Bseq
    #println(cfgnow, SCnow)
    haverec = falses(length(shells))
    wgtsnow = Vector{Float64}(undef, length(shells))
    SCsnow = Vector{ComplexF64}(undef, length(shells))
    sgns = Vector{ComplexF64}(undef, length(shells))
    #每一个观测量是一个元素
    obsers::Vector{Vector{Any}} = []
    #
    resc = real(oldSC)
    #用来记录相对的权重，在一次计算中相对权重一样大
    #
    @push_shell resc shells haverec true
    @push_shell resc shells wgtsnow 1.0
    @push_shell resc shells SCsnow oldSC
    @push_shell resc shells sgns (real(oldSC) / abs(real(oldSC)))
    for (vidx, vexp) in enumerate(valexprs)
        val = oldSC * vexp(Bseq) / abs(real(oldSC))
        push!(obsers, Vector{typeof(val)}(undef, length(shells)))
        @push_shell resc shells obsers[vidx] val
    end
    #
    for itime in Base.OneTo(looptime)
        for ind in CartesianIndices(inicfg)
            prs = @flip oldcfg ind
            prcfg .= oldcfg
            prcfg[ind] = prs
            #
            prwgt = _γ[prs] * oldwgt / _γ[oldcfg[ind]]
            Bseq = Bmarr(prcfg, hub)
            prSC = @sgn prwgt Bseq
            #找到符合其抽样的位置
            isrec = @pull_shell real(prSC) shells haverec
            if isrec
                #如果之前记录过，利用之前记录的进行Metropolis
                SCshell = @pull_shell real(prSC) shells SCsnow
                prpr = min(1.0, abs(real(prSC)) / abs(real(SCshell)))
                if (rand() < prpr)
                    oldwgt = prwgt
                    oldSC = prSC
                    oldcfg .= prcfg
                    #记录到shell
                    resc = real(oldSC)
                    @push_shell resc shells wgtsnow oldwgt
                    @push_shell resc shells SCsnow oldSC 
                end
            else
                #如果之前没有记录，直接加入
                oldwgt = prwgt
                oldSC = prSC
                oldcfg .= prcfg
                #记录到shell
                resc = real(oldSC)
                @push_shell resc shells haverec true
                @push_shell resc shells wgtsnow oldwgt
                @push_shell resc shells SCsnow oldSC
                @push_shell resc shells sgns (real(oldSC) / abs(real(oldSC)))
                for (vidx, vexp) in enumerate(valexprs)
                    val = oldSC * vexp(Bseq) / abs(real(oldSC))
                    @push_shell resc shells obsers[vidx] val
                end
            end
        end
        #进行观测
        @add_shell resc shells sgns (real(oldSC) / abs(real(oldSC)))
        for (vidx, vexp) in enumerate(valexprs)
            val = oldSC * vexp(Bseq) / abs(real(oldSC))
            #if resc < 0
            #    println(oldSC, " ", abs(real(oldSC)), " ", vexp(Bseq), " ", val)
            #end
            @add_shell resc shells obsers[vidx] val
        end
    end
    return sgns, obsers
end


"""
给theta赋值
"""
macro push_theta(sc, shells, nums, dens)
    quote
        __ReSC = real($(esc(sc)))
        for (sidx, sh) in enumerate($(esc(shells)))
            if __ReSC < sh.upper && __ReSC >= sh.lower
                $(esc(nums))[sidx] = __ReSC < sh.middl ? __ReSC / abs(__ReSC) : 0
                $(esc(dens))[sidx] = __ReSC >= sh.middl ? __ReSC / abs(__ReSC) : 0
            end
        end
    end
end


"""
累计
"""
macro add_theta(sc, shells, nums, dens)
    quote
        __ReSC = real($(esc(sc)))
        for (sidx, sh) in enumerate($(esc(shells)))
            if __ReSC < sh.upper && __ReSC >= sh.lower
                $(esc(nums))[sidx] += __ReSC < sh.middl ? __ReSC / abs(__ReSC) : 0
                $(esc(dens))[sidx] += __ReSC >= sh.middl ? __ReSC / abs(__ReSC) : 0
            end
        end
    end
end



"""
omega的流
"""
function omegaflow(looptime::T, 
    hub::GeneralHubbard, inicfg::Matrix{P}, shells::Vector{SignShell}
    ) where {T<:Integer, P<:Integer}
    _γ = [1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    oldcfg = copy(inicfg)
    oldwgt = 1.0
    prcfg = copy(inicfg)
    Bseq = Bmarr(oldcfg, hub)
    oldSC = @sgn 1.0 Bseq
    #println(cfgnow, SCnow)
    haverec = falses(length(shells))
    wgtsnow = zeros(Float64, length(shells))
    SCsnow = zeros(ComplexF64, length(shells))
    numsnow = zeros(ComplexF64,  length(shells))
    densnow = zeros(ComplexF64,  length(shells))
    #
    resc = real(oldSC)
    #用来记录相对的权重，在一次计算中相对权重一样大
    @push_shell resc shells haverec true
    @push_shell resc shells wgtsnow 1.0
    @push_shell resc shells SCsnow oldSC
    #println(resc, numsnow, densnow)
    @push_theta oldSC shells numsnow densnow
    for itime in Base.OneTo(looptime)
        for ind in CartesianIndices(inicfg)
            prs = @flip oldcfg ind
            prcfg .= oldcfg
            prcfg[ind] = prs
            #
            prwgt = _γ[prs] * oldwgt / _γ[oldcfg[ind]]
            Bseq = Bmarr(prcfg, hub)
            prSC = @sgn prwgt Bseq
            #找到符合其抽样的位置
            isrec = @pull_shell real(prSC) shells haverec
            if ismissing(isrec)
                continue
            end
            if isrec
                #如果之前记录过，利用之前记录的进行Metropolis
                SCshell = @pull_shell real(prSC) shells SCsnow
                prpr = min(1.0, abs(real(prSC)) / abs(real(SCshell)))
                if (rand() < prpr)
                    oldwgt = prwgt
                    oldSC = prSC
                    oldcfg .= prcfg
                    #记录到shell
                    resc = real(oldSC)
                    @push_shell resc shells wgtsnow oldwgt
                    @push_shell resc shells SCsnow oldSC 
                end
            else
                #如果之前没有记录，直接加入
                oldwgt = prwgt
                oldSC = prSC
                oldcfg .= prcfg
                #记录到shell
                resc = real(oldSC)
                @push_shell resc shells haverec true
                @push_shell resc shells wgtsnow oldwgt
                @push_shell resc shells SCsnow oldSC
                @push_theta oldSC shells numsnow densnow
            end
        end
        #进行观测
        @add_theta oldSC shells numsnow densnow
    end
    return numsnow, densnow
end


"""
创建观测量的shell
"""
function obser_shells(points::Vector{Float64})
    shells = Vector{SignShell}(undef, length(points)+1)
    shells[1] = SignShell(-Inf, missing, points[1])
    for pidx in Base.OneTo(length(points)-1)
        shells[pidx+1] = SignShell(points[pidx], missing, points[pidx+1])
    end
    shells[end] = SignShell(points[end], missing, Inf)
    return shells
end

"""
创建omega的shell
"""
function omega_shells(shells::Vector{SignShell})
    sh1::Vector{SignShell} = []
    sh2::Vector{SignShell} = []
    olen = length(shells)
    tidx = 1
    while true
        push!(sh1, SignShell(shells[tidx].lower, shells[tidx].upper, shells[tidx+1].upper))
        push!(sh2, SignShell(shells[tidx+1].lower, shells[tidx+1].upper, shells[tidx+2].upper))
        tidx += 2
        if tidx >= olen
            break
        end
    end
    return sh1, sh2
end


"""
计算修正
ob1 = <>^{L0<<Inf}
ob2 = <>^{L1<<L0}
om1 = w(L1<L0<Inf)
om2 = w(L2<L1<L0)
返回
ob3 = <>^{L1<<Inf}
om3 = w(L2<<L1<<Inf)
"""
function omega_step(ob1, om1, ob2, om2)
    ob3 = (om1*ob2 + ob1) / (1 + om1)
    om3 = (om1*om2) / (1 + om1)
    return ob3, om3
end
