#=
利用Monte Carlo进行求和
=#


function mcsum(looptime::T, hub::GeneralHubbard, inicfg::Matrix{P}, 
    valexprs...) where {T<:Integer, P<:Integer}
    _γ = [1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    cfgnow = copy(inicfg)
    prcfg = copy(cfgnow)
    wgtnow = 1.0
    Bseq = Bmarr(cfgnow, hub)
    SCnow = @sgn wgtnow Bseq
    println(cfgnow, SCnow)
    val = []
    sgn = real(SCnow) / abs(real(SCnow))
    for vexp in valexprs
        push!(val, SCnow * vexp(Bseq) / abs(real(SCnow))) 
    end
    #val = zeros(16)#SCnow*wgtnow
    for itime in Base.OneTo(looptime)
        #for (cfg, wgt) in HSIter(2, 1)
        #    Bseq = Bmarr(cfg, $hub)
        #    val += @sgn wgt Bseq
        #end
        for ind in CartesianIndices(inicfg)
            prs = @flip cfgnow ind
            prcfg .= cfgnow
            prcfg[ind] = prs
            #println(ind, " ", prcfg)
            prwgt = _γ[prs] * wgtnow / _γ[cfgnow[ind]]
            Bseq = Bmarr(prcfg, hub)
            prSC = @sgn prwgt Bseq
            prpr = min(1.0, abs(real(prSC)) / abs(real(SCnow)))
            #println(prpr)
            #
            if (rand() < prpr)
                wgtnow = prwgt
                SCnow = prSC
                cfgnow .= prcfg
            end
            #val[(prcfg[1]-1)*4+prcfg[2]] += 1
        end
        sgn += real(SCnow) / abs(real(SCnow))
        for (vidx, vexp) in enumerate(valexprs)
            val[vidx] += SCnow * vexp(Bseq) / abs(real(SCnow))
        end
    end
    return sgn, val
end


