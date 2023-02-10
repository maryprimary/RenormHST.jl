#=
利用Monte Carlo进行求和
=#


macro __mcsum_gfunc(Bseq)
    quote
        siz = size($(esc(Bseq))[1])
        grf = inv(Diagonal(ones(siz[1])) + BprodUDV($(esc(Bseq))))
        grf
    end
end


function mcsum(looptime::T, hub::GeneralHubbard, inicfg::Matrix{P}, 
    valexprs...) where {T<:Integer, P<:Integer}
    _γ = [1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    cfgnow = copy(inicfg)
    wgtnow = 1.0
    Bseqnow = Bmarr(cfgnow, hub)
    Bprodnow = BprodUDV(Bseqnow)
    SCnow = @sgn2 wgtnow Bprodnow
    #创建prcfg的数组
    prcfg = copy(cfgnow)
    println(SCnow)
    #
    val = []
    sgn = real(SCnow) / abs(real(SCnow))
    baregrf = @__mcsum_gfunc Bseqnow
    grf = SCnow / (abs(real(SCnow))) * baregrf
    for vexp in valexprs
        push!(val, SCnow * vexp(baregrf) / abs(real(SCnow))) 
    end
    #val = zeros(16)#SCnow*wgtnow
    for itime in Base.OneTo(looptime)
        #for (cfg, wgt) in HSIter(2, 1)
        #    Bseq = Bmarr(cfg, $hub)
        #    val += @sgn wgt Bseq
        #end
        # 重置参数
        nowtau = 1
        if mod(itime, update_control.Bpdre) == 0
            Bprodnow = BprodUDV(Bseqnow)
        end
        for ind in CartesianIndices(inicfg)
            #如果此时已经超过了, 推进Bseq和Bprod
            if nowtau != ind[1]
                nowtau = ind[1]
                #println(Bprodnow)
                #println(BprodUDV(Bseqnow))
                #推进Bprod, 要优先推进Bprod，因为这里用到了Bseqnow[1]
                #TODO: 将Bprod和Bseq放在一起推进
                Bprodnow = @uptau_Bprod! Bprodnow Bseqnow
                #推进Bseq
                Bseqnow = @uptau_Bseq! Bseqnow
            end
            prs = @flip cfgnow ind
            prcfg .= cfgnow
            prcfg[ind] = prs
            #println(ind, " ", prcfg)
            prwgt = _γ[prs] * wgtnow / _γ[cfgnow[ind]]
            #使用Bmarr
            #prBseq = Bmarr(prcfg, hub)
            #使用迭代Btop
            #prBseq = @hspr_Bseq Bseqnow hub prcfg[nowtau, :]
            prBseq = @hspr_Bmat Bseqnow hub ind[2] prcfg[ind] cfgnow[ind]
            #迭代Bprod
            prBprod = @hspr_Bprod Bprodnow hub ind[2] prcfg[ind] cfgnow[ind]
            #
            #prSC = @sgn prwgt prBseq
            prSC = @sgn2 prwgt prBprod
            prpr = min(1.0, abs(real(prSC)) / abs(real(SCnow)))
            #println(prpr)
            #
            if (rand() < prpr)
                wgtnow = prwgt
                SCnow = prSC
                cfgnow .= prcfg
                Bseqnow = prBseq
                Bprodnow = prBprod
            end
        end
        #重新把1放到前面, 先Bprod，后Bseq
        @uptau_Bprod! Bprodnow Bseqnow
        @uptau_Bseq! Bseqnow
        #使用Bmarr
        #Bseqnow = Bmarr(cfgnow, hub)
        #
        sgn += real(SCnow) / abs(real(SCnow))
        baregrf = @__mcsum_gfunc Bseqnow
        grf += SCnow * baregrf / abs(real(SCnow))
        for (vidx, vexp) in enumerate(valexprs)
            val[vidx] += SCnow * vexp(baregrf) / abs(real(SCnow))
        end
    end
    return sgn, grf, val
end


