function    beta=get_beta(Q,Aeq,Aiq,mystruct,Qc)
    Q=FindQR(Q,Aeq,Aiq,mystruct);
    beta=linsolve(Q,Qc);
end