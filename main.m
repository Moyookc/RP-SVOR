function [out]=main(data_flag)
    global KTYPE KSCALE fake_zero size_training fake_zero_linprog NumSingular
    % KTYPE = 6;
    fake_zero=10^-8;
    fake_zero2=10^-8;
    NumSingular=0;
    t_range=[0,10000000];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [xx,yy,Q,c,q,Aeq,beq,peq,Aiq,biq,piq,L,U,O]=PrepareMatsORwei(t_range,data_flag);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out=GSP(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,L,U,O,t_range,xx,yy);
end
