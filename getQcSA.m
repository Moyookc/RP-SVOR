function Qc=getQcSA(mystruct,Q,q,Aeq,peq,Aiq,piq,O)
    SE=mystruct.SE;
    SR=mystruct.SR;
    indexJ=mystruct.indexJ;
    p=[piq;peq];
    AA=[Aiq Aeq]';
    Qc=[-q(SE)-Q(SE,SR)*O(SR)];
    Qc=[Qc;p(indexJ)-AA(indexJ,SR)*O(SR)];
end