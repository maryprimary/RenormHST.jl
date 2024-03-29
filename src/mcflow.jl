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


macro __mcflow_gfunc(Bseq)
    quote
        siz = size($(esc(Bseq))[1])
        grf = inv(Diagonal(ones(siz[1])) + BprodUDV($(esc(Bseq))))
        grf
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
    oldBseq = Bmarr(oldcfg, hub)
    #oldSC = @sgn 1.0 oldBseq
    oldBprod = BprodUDV(oldBseq)
    oldSC = @sgn2 1.0 oldBprod
    #@sgn 1.0 oldBseq
    #创建prcfg的数组
    prcfg = copy(inicfg)
    #println(cfgnow, SCnow)
    haverec = falses(length(shells))
    wgtsnow = Vector{Float64}(undef, length(shells))
    SCsnow = Vector{ComplexF64}(undef, length(shells))
    sgns = Vector{ComplexF64}(undef, length(shells))
    grfs = Vector{Matrix{ComplexF64}}(undef, length(shells))
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
    baregrf = @__mcflow_gfunc oldBseq
    @push_shell resc shells grfs oldSC * baregrf / (abs(real(oldSC)))
    for (vidx, vexp) in enumerate(valexprs)
        val = oldSC * vexp(baregrf) / abs(real(oldSC))
        push!(obsers, Vector{typeof(val)}(undef, length(shells)))
        @push_shell resc shells obsers[vidx] val
    end
    #
    for itime in Base.OneTo(looptime)
        #重置nowtau
        nowtau = 1
        if mod(itime, update_control.Bpdre) == 0
            oldBprod = BprodUDV(oldBseq)
        end
        for ind in CartesianIndices(inicfg)
            if nowtau != ind[1]
                nowtau = ind[1]
                #推进Bprod
                oldBprod = @uptau_Bprod! oldBprod oldBseq
                #推进Bseq
                oldBseq = @uptau_Bseq! oldBseq
            end
            prs = @flip oldcfg ind
            prcfg .= oldcfg
            prcfg[ind] = prs
            #
            prwgt = _γ[prs] * oldwgt / _γ[oldcfg[ind]]
            #使用Bmarr
            #prBseq = Bmarr(prcfg, hub)
            #prSC = @sgn prwgt prBseq
            #使用prBseq
            prBseq = @hspr_Bmat oldBseq hub ind[2] prcfg[ind] oldcfg[ind]
            #prSC = @sgn prwgt prBseq
            prBprod = @hspr_Bprod oldBprod hub ind[2] prcfg[ind] oldcfg[ind]
            prSC = @sgn2 prwgt prBprod
            #找到符合其抽样的位置
            #isrec = @pull_shell real(prSC) shells haverec
            #如果查找不到，代表没有shell能包含现在的HS
            #直接continue可能导致一直没有合适的HS
            #if ismissing(isrec)
            #    #println(real(prSC), " ", shells)
            #    continue
            #end
            #if isrec
                #如果之前记录过，利用之前记录的进行Metropolis
                #SCshell = @pull_shell real(prSC) shells SCsnow
            prpr = min(1.0, abs(real(prSC)) / abs(real(oldSC)))
            if (rand() < prpr)
                oldwgt = prwgt
                oldSC = prSC
                oldcfg .= prcfg
                oldBseq = prBseq
                oldBprod = prBprod
                #记录到shell
                resc = real(oldSC)
                @push_shell resc shells wgtsnow oldwgt
                @push_shell resc shells SCsnow oldSC 
            end
        end
        #重新把1放到前面
        @uptau_Bprod! oldBprod oldBseq
        @uptau_Bseq! oldBseq
        #使用Bmarr
        #oldBseq = Bmarr(oldcfg, hub)
        #println(oldBseq)
        #throw(error())
        resc = real(oldSC)
        isrec = @pull_shell resc shells haverec
        if ismissing(isrec)
            continue
        end
        #进行观测
        if isrec
            @add_shell resc shells sgns (real(oldSC) / abs(real(oldSC)))
            baregrf = @__mcflow_gfunc oldBseq
            @add_shell resc shells grfs oldSC * baregrf / abs(real(oldSC))
            for (vidx, vexp) in enumerate(valexprs)
                val = oldSC * vexp(baregrf) / abs(real(oldSC))
                #if resc < 0
                #    println(oldSC, " ", abs(real(oldSC)), " ", vexp(Bseq), " ", val)
                #end
                @add_shell resc shells obsers[vidx] val
            end
        else
            @push_shell resc shells haverec true
            @push_shell resc shells sgns (real(oldSC) / abs(real(oldSC)))
            baregrf = @__mcflow_gfunc oldBseq
            @push_shell resc shells grfs oldSC * baregrf / abs(real(oldSC))
            for (vidx, vexp) in enumerate(valexprs)
                val = oldSC * vexp(baregrf) / abs(real(oldSC))
                #if resc < 0
                #    println(oldSC, " ", abs(real(oldSC)), " ", vexp(Bseq), " ", val)
                #end
                @push_shell resc shells obsers[vidx] val
            end
        end
    end
    return sgns, grfs, obsers
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
    oldBseq = Bmarr(oldcfg, hub)
    #oldBprod = BprodUDV(oldBseq)
    oldSC = @sgn2 1.0 oldBprod
    #创建prcfg的数组
    prcfg = copy(inicfg)
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
    #prcount = 0
    #prpos = 0
    for itime in Base.OneTo(looptime)
        #重置nowtau
        nowtau = 1
        if mod(itime, update_control.Bpdre) == 0
            #oldBprod = BprodUDV(oldBseq)
        end
        for ind in CartesianIndices(inicfg)
            if nowtau != ind[1]
                nowtau = ind[1]
                #推进Bprod
                #oldBprod = @uptau_Bprod! oldBprod oldBseq
                #推进Bseq
                oldBseq = @uptau_Bseq! oldBseq
            end
            prs = @flip oldcfg ind
            prcfg .= oldcfg
            prcfg[ind] = prs
            #
            prwgt = _γ[prs] * oldwgt / _γ[oldcfg[ind]]
            #使用Bmarr
            #prBseq = Bmarr(prcfg, hub)
            #prSC = @sgn prwgt prBseq
            #迭代
            prBseq = @hspr_Bmat oldBseq hub ind[2] prcfg[ind] oldcfg[ind]
            #prBprod = @hspr_Bprod oldBprod hub ind[2] prcfg[ind] oldcfg[ind]
            prSC = @sgn2 prwgt prBprod
            #找到符合其抽样的位置
            #isrec = @pull_shell real(prSC) shells haverec
            #println(isrec)
            #如果查找不到，代表没有shell能包含现在的HS
            #直接continue可能导致一直没有合适的HS
            #如果之前记录过，利用之前记录的进行Metropolis
            prpr = min(1.0, abs(real(prSC)) / abs(real(oldSC)))
            if (rand() < prpr)
                oldwgt = prwgt
                oldSC = prSC
                oldcfg .= prcfg
                oldBseq = prBseq
                #oldBprod = prBprod
                #记录到shell
                resc = real(oldSC)
                @push_shell resc shells wgtsnow oldwgt
                @push_shell resc shells SCsnow oldSC 
                #prcount += 1
                #println(oldwgt)
            end
        end
        #重新把tau变成1
        @uptau_Bprod! oldBprod oldBseq
        @uptau_Bseq! oldBseq
        #bseq2 = Bmarr(oldcfg, hub)
        #println((@sgn 1.0 bseq2), "bmarr")
        #println((@sgn2 1.0 oldBprod), "bprod")
        #进行观测
        resc = real(oldSC)
        isrec = @pull_shell resc shells haverec
        if ismissing(isrec)
            continue
        end
        if isrec
            @add_theta oldSC shells numsnow densnow
        else
            @push_shell resc shells haverec true
            @push_theta oldSC shells numsnow densnow
        end
    end
    #println("omega flow, prcount ", prcount, " ", shells, " ", prpos) 
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
        if (tidx+2 <= length(shells))
            push!(sh2, SignShell(shells[tidx+1].lower, shells[tidx+1].upper, shells[tidx+2].upper))
        end
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


"""
计算omega所有的
"""
function omega_steps(oshell, ompris, omeprs)
    mshell = Vector{SignShell}(undef, length(oshell))
    omegal = zeros(ComplexF64, length(oshell))
    omerrl = zeros(ComplexF64, length(oshell))
    mshell[end] = SignShell(oshell[end].lower, oshell[end].middl, oshell[end].upper) 
    omegal[end] = ompris[end]
    omerrl[end] = omeprs[end]
    for oidx in Base.OneTo(length(oshell)-1)
        if (mshell[end-oidx+1].lower != oshell[end-oidx].middl ||
            mshell[end-oidx+1].middl != oshell[end-oidx].upper)
            throw(error("oshell "*string(oidx)*"do not match"))
        end
        mshell[end-oidx] = SignShell(oshell[end-oidx].lower, oshell[end-oidx].middl,
        mshell[end-oidx+1].upper)
        omegal[end-oidx] = ompris[end-oidx]*omegal[end-oidx+1] / (1 + omegal[end-oidx+1])
        errnum2 = ^(ompris[end-oidx], 2) * ^(omerrl[end-oidx+1], 2) +
        ^(omeprs[end-oidx], 2) * ^(omegal[end-oidx+1], 2)
        omerrl[end-oidx] = errnum2 / ^(1+omegal[end-oidx+1], 2) +
        ^(omerrl[end-oidx+1], 2) * ^(ompris[end-oidx]*omegal[end-oidx+1], 2) / ^(1+omegal[end-oidx+1], 4)
        omerrl[end-oidx] = sqrt(omerrl[end-oidx]) 
    end
    #mshell[1] = SignShell(oshell[1].lower, oshell[1].middl, mshell[2].upper)
    #omegal[1] = ompris[1]*omegal[2] / (1 + omegal[2])
    return mshell, omegal, omerrl
end


"""
计算所有的obser
"""
function obser_steps(shells, obsers, momegas, obser_err, momega_err)
    mobser = Vector{typeof(obsers[end])}(undef, length(shells))
    mobser_err = Vector{typeof(obsers[end])}(undef, length(shells))
    mshell2 = Vector{SignShell}(undef, length(shells))
    mobser[end] = obsers[end]
    mobser_err[end] = obser_err[end]
    mshell2[end] = SignShell(shells[end].lower, missing, shells[end].upper)
    for oidx in Base.OneTo(length(shells)-1)
        if (shells[end-oidx].upper != mshell2[end-oidx+1].lower)
            throw(error("shells "*string(oidx)*"do not match"))
        end
        mshell2[end-oidx] = SignShell(shells[end-oidx].lower, missing, mshell2[end-oidx+1].upper)
        mobser[end-oidx] = (momegas[end-oidx+1]*obsers[end-oidx] + mobser[end-oidx+1]) / (1 + momegas[end-oidx+1])
        err1 = ^(momega_err[end-oidx+1], 2)*(.^(obsers[end-oidx], 2)) +
        ^(momegas[end-oidx+1], 2)*(.^(obser_err[end-oidx], 2))
        println(oidx, " 1:", err1)
        err1 = sqrt(err1) + mobser_err[end-oidx+1] 
        println(oidx, " 2:", err1)
        err1 = .^(err1, 2) / ^(1+momegas[end-oidx+1], 2) +
        ^(momega_err[end-oidx+1], 2)*(.^(momegas[end-oidx+1]*obsers[end-oidx] + mobser[end-oidx+1], 2)) /
        ^(1+momegas[end-oidx+1], 4)
        println(oidx, " 3:", err1)
        mobser_err[end-oidx] = sqrt(err1)
    end
    return mshell2, mobser, mobser_err
end


