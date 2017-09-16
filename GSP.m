function out=GSP(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,L,U,O,t_range,x,y)
global KTYPE KSCALE DeltaNu fake_zero 
    t=t_range(1);
    tic
    path=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%core code
        initial_solution=initial(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,L,U,O,t_range,x,y);
        [local_path]=SolutionPath(initial_solution,Q,c,q,Aeq,beq,peq,Aiq,biq,piq,L,U,O,t_range);
        path=[path;local_path];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc
    time=toc ;
    steps=length(path);
    out.Steps=steps;
    out.Path=path;
end