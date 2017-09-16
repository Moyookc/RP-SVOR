function [max_change_alpha index_max_change_alpha]=...
    find_positive_change_alphaSA(mystruct,beta,L,U,O)
    global fake_zero 
    t=mystruct.t;
    max_change_alpha=inf;
    index_max_change_alpha=inf;
    support_set=mystruct.SE;
    alpha = mystruct.alpha;
    local_alpha=alpha(support_set);
    len1=length(support_set);
    local_beta=beta([1:len1]);
    local_O=O(support_set);
    local_index=find((local_beta<fake_zero & local_beta>-fake_zero) | (local_beta-local_O)<0);
    local_alpha(local_index)=[];
    support_set(local_index)=[];
    local_beta(local_index)=[];
    if length(support_set)>0
        local_O=O(support_set);
        local_C=(U(support_set)+local_O*t).*(local_beta>0)+(L(support_set)-local_alpha).*(local_beta<0);
        local_bb=local_beta.*(local_beta<0)+(local_beta-local_O).*(local_beta>0)
        local_change_al=local_C./local_bb;
        [min_alpha_change min_alpha_change_index]=min(local_change_al);
        if length(local_change_al)==0
            max_change_alpha=inf;
        else
            max_change_alpha=min_alpha_change;
            index_max_change_alpha=support_set(min_alpha_change_index);
        end
    end
end