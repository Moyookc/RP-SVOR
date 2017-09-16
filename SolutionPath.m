function [path]=SolutionPath(initial_solution,Q,c,q,Aeq,beq,peq,Aiq,biq,piq,L,U,O,t_range)
    global fake_zero
    path=[initial_solution];
    mystruct=initial_solution;
    Qc=getQcSA(mystruct,Q,q,Aeq,peq,Aiq,piq,O);
    [g, h]=get_g(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,mystruct);
    steps=0;
    while mystruct.t<(t_range(2)-fake_zero)
        beta=get_beta(Q,Aeq,Aiq,mystruct,Qc);
        [gama_g, gama_h]=get_gamaSA(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,O,mystruct,beta);
        [max_delta_zeta,index]=get_max_changeSA(mystruct,g, h,beta,gama_g, gama_h,L, U, O, Aiq, t_range);
        [g,h,mystruct,Qc]=updata(mystruct,g,h,beta,gama_g, gama_h,max_delta_zeta,index,Q,q,Aeq,peq,Aiq,piq,L,U,O);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%accessorial code
        out_KKT1=test_KKT(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,mystruct);
        mystruct.KKT=out_KKT1;
        path=[path;mystruct];
        steps=steps+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end