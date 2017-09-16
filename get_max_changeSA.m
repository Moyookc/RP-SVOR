function [max_delta_zeta,index]=get_max_changeSA(mystruct,g, h,beta,gama_g, gama_h, L, U, O, Aiq, t_range)
    [max_change_alpha index_max_change_alpha]=find_positive_change_alphaSA...
        (mystruct,beta,L,U,O);
    [max_change_u index_max_change_u]=find_positive_change_u...
        (mystruct,beta,L,U,O,Aiq);
    [max_change_g index_max_change_g]=find_max_change_g...
        (mystruct,g,gama_g);
    [max_change_h index_max_change_h]=find_max_change_h...
        (mystruct,h,gama_h);
    [max_change_sum_alpha]=find_change_sum_alpha...
        (mystruct,t_range);
    [max_delta_zeta index]=find_max_change...
        (max_change_alpha, index_max_change_alpha,max_change_u,index_max_change_u,max_change_g,...
        index_max_change_g,max_change_h,index_max_change_h,max_change_sum_alpha);
end
function   [max_change_u index_max_change_u]=find_positive_change_u...
        (mystruct,beta,L,U,O,Aiq)
    global fake_zero
    t=mystruct.t;
    support_set=mystruct.SE;
    len1=length(support_set);
    indexJ=mystruct.indexJ;
    if length(Aiq)==0
        tmp_len2=0;
    else
        tmp_len2=length(Aiq(1,:));
    end  
    max_change_u=inf;
    index_max_change_u=inf;
    local_indexJ=indexJ;
    lang_u = mystruct.lang_u;
    [a,b]=find(indexJ>(tmp_len2));
    indexJ(a)=[];
    if length(indexJ)>0
        local_beta_u=beta([len1+1:end]);
        local_beta_u=local_beta_u(indexJ);
        local_lang_u=lang_u(indexJ);
        local_index=find(local_beta_u>-fake_zero);
        local_beta_u(local_index)=[];
        local_indexJ(local_index)=[];
        local_lang_u(local_index)=[];
        local_change_al=-local_lang_u./local_beta_u;
        [min_alpha_change min_alpha_change_index]=min(local_change_al);
        if length(local_change_al)==0
            max_change_alpha=inf;
        else
            max_change_alpha=min_alpha_change;
            index_max_change_alpha=local_indexJ(min_alpha_change_index);
        end
    end
end
function [max_change_h index_max_change_h]=find_max_change_h(mystruct,h,gama)
    global fake_zero
    other_set=mystruct.indexCJ;
    gama=gama(other_set);
    max_change_h=inf;
    index_max_change_h=inf;
    local_index=find(gama<fake_zero & gama>-fake_zero);
    gama(local_index)=[];
    other_set(local_index)=[];
    local_h=h(other_set);
    local_index2=find(local_h<fake_zero & local_h>-fake_zero);
    local_h(local_index2)=[];
    gama(local_index2)=[];
    other_set(local_index2)=[];
    local_index3=find(local_h.*gama>=0);
    local_h(local_index3)=[];
    gama(local_index3)=[];
    other_set(local_index3)=[];
    local_change_h=-local_h./gama;
    [local_change local_index]=min(local_change_h);
    if length(local_change_h)==0
        max_change_h=inf;
    else
        max_change_h=local_change;
        index_max_change_h=other_set(local_index);
    end
end