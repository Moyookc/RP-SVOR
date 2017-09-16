function out_KKT=sub_KKT_g(g,mystruct)
    global fake_zero
    out_KKT=1;
    SE=mystruct.SE;
    SR=mystruct.SR;
    SL=mystruct.SL;
    tmp1=sum(g(SE)>fake_zero)+sum(g(SE)<-fake_zero);
    if tmp1>0
        out_KKT=0;
        return;
    end
    tmp2=sum(g(SR)>fake_zero);
    if tmp2>0
        out_KKT=0;
        return;
    end
    tmp3=sum(g(SL)<-fake_zero);
    if tmp3>0
        out_KKT=0;
        return;
    end
end