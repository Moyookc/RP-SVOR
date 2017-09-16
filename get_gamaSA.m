function    [gama_g, gama_h]=get_gamaSA(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,O,mystruct,beta)
    SE=mystruct.SE;
    SR=mystruct.SR;
    indexJ=mystruct.indexJ;
    AA=[Aiq Aeq];
    pp=[piq;peq];
    len_SE=length(SE);
    alpha=mystruct.alpha;
    beta_alpha=beta([1:len_SE]);
    beta_u=beta([len_SE+1:end]);
    gama_g = Q(:,SE)*beta_alpha+Q(:,SR)*O(SR)+q+AA(:,indexJ)*beta_u;
    gama_h = AA(SE,:)'*beta_alpha+AA(SR,:)'*O(SR) - pp;
end