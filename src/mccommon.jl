#=
mc通用的一些方法
=#


function warmup_configuration(looptime::T, hub::GeneralHubbard, inicfg::Matrix{P}
    ) where {T<:Integer, P<:Integer}
    _γ = [1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    cfgnow = copy(inicfg)
    wgtnow = 1.0
    Bseqnow = Bmarr(cfgnow, hub)
    Bprodnow = BprodUDV(Bseqnow)
    SCnow = @sgn2 wgtnow Bprodnow
    #创建prcfg的数组
    prcfg = copy(inicfg)
    #
    for itime in Base.OneTo(looptime)
        nowtau = 1
        if mod(itime, update_control.Bpdre) == 0
            Bprodnow = BprodUDV(Bseqnow)
        end
        for ind in CartesianIndices(inicfg)
            #如果此时已经超过了, 推进Bseq和Bprod
            if nowtau != ind[1]
                nowtau = ind[1]
                #推进Bprod, 要优先推进Bprod，因为这里用到了Bseqnow[1]
                #TODO: 将Bprod和Bseq放在一起推进
                Bprodnow = @uptau_Bprod! Bprodnow Bseqnow
                #推进Bseq
                Bseqnow = @uptau_Bseq! Bseqnow
            end
            prs = @flip cfgnow ind
            prcfg .= cfgnow
            prcfg[ind] = prs
            #
            prwgt = _γ[prs] * wgtnow / _γ[cfgnow[ind]]
            #迭代Bseq和Bprod
            prBseq = @hspr_Bmat Bseqnow hub ind[2] prcfg[ind] cfgnow[ind]
            prBprod = @hspr_Bprod Bprodnow hub ind[2] prcfg[ind] cfgnow[ind]
            prSC = @sgn2 prwgt prBprod
            prpr = min(1.0, abs(real(prSC)) / abs(real(SCnow)))
            if (rand() < prpr)
                wgtnow = prwgt
                SCnow = prSC
                cfgnow .= prcfg
                Bseqnow = prBseq
                Bprodnow = prBprod
            end
        end
        #重新把1放到前面
        @uptau_Bprod! Bprodnow Bseqnow
        @uptau_Bseq! Bseqnow
    end
    return cfgnow
end



"""
统计符号的情况
"""
function sgn_hist(looptime::T, hub::GeneralHubbard, inicfg::Matrix{P}, maxsgn::Int64
    ) where {T<:Integer, P<:Integer}
    _γ = [1-(√6/3), 1+(√6/3), 1+(√6/3), 1-(√6/3)]
    cfgnow = copy(inicfg)
    wgtnow = 1.0
    Bseqnow = Bmarr(cfgnow, hub)
    Bprodnow = BprodUDV(Bseqnow)
    SCnow = @sgn2 wgtnow Bprodnow
    #
    prcfg = copy(inicfg)
    hist = zeros(Int64, maxsgn+1)
    #
    for itime in Base.OneTo(looptime)
        nowtau = 1
        if mod(itime, update_control.Bpdre) == 0
            Bprodnow = BprodUDV(Bseqnow)
        end
        for ind in CartesianIndices(inicfg)
            #如果此时已经超过了, 推进Bseq和Bprod
            if nowtau != ind[1]
                nowtau = ind[1]
                #推进Bprod, 要优先推进Bprod，因为这里用到了Bseqnow[1]
                #TODO: 将Bprod和Bseq放在一起推进
                Bprodnow = @uptau_Bprod! Bprodnow Bseqnow
                #推进Bseq
                Bseqnow = @uptau_Bseq! Bseqnow
            end
            prcfg .= cfgnow
            oldprs = cfgnow[ind]
            prs = @flip cfgnow ind
            prcfg[ind] = prs
            #
            prwgt = _γ[prs] * wgtnow / _γ[oldprs]
            #迭代Bseq和Bprod
            prBseq = @hspr_Bmat Bseqnow hub ind[2] prcfg[ind] oldprs
            prBprod = @hspr_Bprod Bprodnow hub ind[2] prcfg[ind] oldprs
            prSC = @sgn2 prwgt prBprod
            prpr = min(1.0, abs(real(prSC)) / abs(real(SCnow)))
            if (rand() < prpr)
                wgtnow = prwgt
                SCnow = prSC
                cfgnow .= prcfg
                Bseqnow = prBseq
                Bprodnow = prBprod            
            end
        end
        #重新把1放到前面
        @uptau_Bprod! Bprodnow Bseqnow
        @uptau_Bseq! Bseqnow
        #观测
        resc = real(SCnow)
        ridx = min(0.0, resc)
        ridx = min(maxsgn, -ridx)
        ridx = Int64(ceil(ridx))
        hist[ridx+1] += 1
    end
    return hist
end


