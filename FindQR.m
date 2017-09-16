function Q=FindQR(Q,Aeq,Aiq,mystruct)
global NumSingular fake_zero
    SE=mystruct.SE;
    SR=mystruct.SR;
    indexJ=mystruct.indexJ;
    len1=length(indexJ);
    AA=[Aiq Aeq]';
    Q=[Q(SE,SE) AA(indexJ,SE)';AA(indexJ,SE) zeros(len1,len1)];
    tmp=det(Q);
    if tmp<fake_zero & tmp>-fake_zero
        NumSingular=NumSingular+1;
    end
end