function out_KKT=test_KKT(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,mystruct)
    out_KKT=1;
    [g, h]=get_g(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,mystruct);
    local_tmp=sub_KKT_g(g,mystruct);
    if local_tmp==0
        out_KKT=0;
        return;
    end
    local_tmp=sub_KKT_h(h,mystruct);
    if local_tmp==0
        out_KKT=0;
        return;
    end
end
function out_KKT=sub_KKT_h(h,mystruct)
    global fake_zero
    out_KKT=1;
    indexJ=mystruct.indexJ;
    indexCJ=mystruct.indexCJ;
    tmp1=sum(h(indexJ)>fake_zero)+sum(h(indexJ)<-fake_zero);
    if tmp1>0
        out_KKT=0;
        return;
    end
    tmp2=sum(h(indexCJ)>fake_zero);
    if tmp2>0
        out_KKT=0;
        return;
    end
end