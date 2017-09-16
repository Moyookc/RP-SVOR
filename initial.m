function out=initial(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,L,U,O,t_range,xx,yy)
    global solver_type
    t=t_range(1);
    cc=c+q*t;
    bbiq=biq+piq*t;
    bbeq=beq+peq*t;
    UU=U+O*t;
    if solver_type==1
    %     options = optimoptions('quadprog',...
    %     'Algorithm','interior-point-convex','TolX',1e-16,'TolFun',1e-16,'MaxIter',10000);
        options = optimoptions('quadprog',...
        'Algorithm','active-set','TolX',1e-20,'TolFun',1e-20,'MaxIter',2000);
        [alpha,fval,exitflag,output]=quadprog(Q,cc,Aiq',bbiq,Aeq',bbeq,L,UU,[],options);
        if exitflag==1
            [mystruct]=GetInitialStruct(alpha,Q,cc,Aiq,bbiq,Aeq,bbeq,L,UU,t);
            out=mystruct;
            out_KKT=test_KKT(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,mystruct);
            out.KKT=out_KKT;
        else
            error('quadprog cannot return a optimal solution!!');
        end
    end
    if solver_type==2
        C=UU(1);
        out_svor=SVOR(xx,yy,C);
        mystruct=GetInitialStructSmo(out_svor,L,UU,t);
        out=mystruct;
        out_KKT=test_KKT(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,mystruct);
        out.KKT=out_KKT;
    end
end
function mystruct=GetInitialStructSmo(out_svor,L,UU,t);
    global fake_zero size_training
    mystruct=[];
    size_subclass=out_svor.size_subclass;
    alpha=out_svor.alpha;
    length_classes=length(out_svor.b)+1;
    local_index=[];
    index1=0;
    for i=1:length_classes-1
        local_index=[local_index;[index1+size_subclass(i)+1:index1+size_subclass(i)+size_subclass(i+1)]'];
        index1=index1+size_subclass(i)+size_subclass(i+1);
    end
    local_apha=alpha(local_index);
    alpha(local_index)=[];
    alpha=[alpha;local_apha];
    alpha=[alpha;out_svor.mu];
    SE=find(alpha>(L+fake_zero) &  alpha<(UU-fake_zero));
    SL=find(alpha<=(L+fake_zero));
    SR=find(alpha>=(UU-fake_zero));
    indexJ=[1:length_classes-1];
    indexCJ=[];
    mystruct.alpha=alpha;
    mystruct.t=t;
    mystruct.SE=SE;
    mystruct.SL=SL;
    mystruct.SR=SR;
    mystruct.indexJ=indexJ;
    mystruct.indexCJ=indexCJ;
    mystruct.lang_u=out_svor.b;
end
