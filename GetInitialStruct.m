function [mystruct]=GetInitialStruct(alpha,Q,cc,Aiq,bbiq,Aeq,bbeq,L,UU,t)
    global fake_zero size_training
    mystruct=[];
    mystruct.alpha=alpha;
    mystruct.t=t;
    [lang_u,SE,SL,SR,indexJ,indexCJ]=GetStructElements(alpha,Q,cc,Aiq,bbiq,Aeq,bbeq,L,UU);
    mystruct.SE=SE;
    mystruct.SL=SL;
    mystruct.SR=SR;
    mystruct.indexJ=indexJ;
    mystruct.indexCJ=indexCJ;
    mystruct.lang_u=lang_u;
end

