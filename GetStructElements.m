function [lang_u,SE,SL,SR,indexJ,indexCJ]=GetStructElements(alpha,Q,cc,Aiq,bbiq,Aeq,bbeq,L,UU)
    global fake_zero fake_zero_linprog
    fake_zero_linprog=10^-8;
    SE=find(alpha>(L+fake_zero) &  alpha<(UU-fake_zero));
    SL=find(alpha<=(L+fake_zero));
    SR=find(alpha>=(UU-fake_zero));
    tmp1=Q*alpha+cc;
    AA=[Aiq';Aeq'];
    options = optimset('Algorithm','active-set','MaxIter',500,'TolFun',1e-12);
    if length(SE)<2
        lin_Aeq=AA(:,SE)';
        lin_beq=-tmp1(SE);
        mylen=length(AA(:,1));
        lin_f=ones(mylen,1);
        if length(Aiq)>0
            tmp_len1=length(Aiq(1,:));
        else
            tmp_len1=0;
        end
        if length(Aeq)>0
            tmp_len2=length(Aeq(1,:));
        else
            tmp_len2=0;
        end
        lin_lb=[zeros(tmp_len1,1);-inf*ones(tmp_len2,1)];
        lin_ub=inf*ones(mylen,1);
        lin_A=[-AA(:,SL) AA(:,SR)]';
        lin_b=[tmp1(SL);-tmp1(SR)];
        if length(lin_A)==0
            lin_A=[];
        end
        if length(lin_b)==0
            lin_b=[];
        end
        [lang_u, fval, exitflag,output]= linprog(lin_f,lin_A,lin_b,lin_Aeq,lin_beq,lin_lb,lin_ub,[],options);
    else
        mylen=length(AA(:,1));
        lin_f=ones(mylen,1);
        if length(Aiq)>0
            tmp_len1=length(Aiq(1,:));
        else
            tmp_len1=0;
        end
        if length(Aeq)>0
            tmp_len2=length(Aeq(1,:));
        else
            tmp_len2=0;
        end
        lin_lb=[zeros(tmp_len1,1);-inf*ones(tmp_len2,1)];
        lin_ub=inf*ones(mylen,1);
        lin_A=[-AA(:,SL) -AA(:,SE) AA(:,SR) AA(:,SE)]';
        lin_b=[tmp1(SL);tmp1(SE)+fake_zero_linprog;-tmp1(SR);-tmp1(SE)+fake_zero_linprog];
        if length(lin_A)==0
            lin_A=[];
        end
        if length(lin_b)==0
            lin_b=[];
        end
        lin_Aeq=[];
        lin_beq=[];
        [lang_u, fval, exitflag,output]= linprog(lin_f,lin_A,lin_b,lin_Aeq,lin_beq,lin_lb,lin_ub,[],options);
    end
    if exitflag==1
        if length(Aiq)==0
            index_iq=[];
        else
            index_iq=[1:length(Aiq(1,:))];
        end
        if length(Aeq)==0
            index_eq=[];
        else
            index_eq=[length(index_iq)+1:length(index_iq)+length(Aeq(1,:))]';
        end
        lang_u_iq=lang_u(index_iq);
        if length(lang_u_iq)==0
            indexJ=index_eq;
            indexCJ=[];
        else
            indexJ=[find(lang_u_iq >= fake_zero);index_eq];
            indexCJ=find(lang_u_iq < fake_zero);
        end
    else
        error('linprog cannot find a feasible solution!!!')
    end
end