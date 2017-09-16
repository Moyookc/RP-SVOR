function [max_change_g index_max_change_g]=find_max_change_g(mystruct,g,gama)
    global fake_zero
    SL=mystruct.SL;
    SR=mystruct.SR;
    other_set=sort([SL;SR]);
    gama=gama(other_set);
    max_change_g=inf;
    index_max_change_g=inf;
    local_index=find(gama<fake_zero & gama>-fake_zero);
    gama(local_index)=[];
    other_set(local_index)=[];
    local_g=g(other_set);
    local_index2=find(local_g<fake_zero & local_g>-fake_zero);
    local_g(local_index2)=[];
    gama(local_index2)=[];
    other_set(local_index2)=[];
    local_index3=find(local_g.*gama>=0);
    local_g(local_index3)=[];
    gama(local_index3)=[];
    other_set(local_index3)=[];
    local_change_g=-local_g./gama;
    [local_change local_index]=min(local_change_g);
    if length(local_change_g)==0
        max_change_g=inf;
    else
        max_change_g=local_change;
        index_max_change_g=other_set(local_index);
    end
end