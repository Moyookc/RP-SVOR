function [g, h]=get_g(Q,c,q,Aeq,beq,peq,Aiq,biq,piq,mystruct)
    t=mystruct.t;
    alpha=mystruct.alpha;
    lang_u=mystruct.lang_u;
    AA=[Aiq Aeq]';
    bb=[biq;beq];
    pp=[piq;peq];
    g=Q*alpha+c+q*t+AA'*lang_u;
    h=AA*alpha-bb-pp*t;
end