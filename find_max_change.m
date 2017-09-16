function [updata_max_change updata_index]=find_max_change...
    (max_change_alpha, index_max_change_alpha,max_change_u,index_max_change_u,max_change_g,...
        index_max_change_g,max_change_h,index_max_change_h,max_change_sum_alpha)
    updata_max_change=inf;
    updata_index=0;
    if updata_max_change>max_change_alpha
        updata_max_change=max_change_alpha;
        updata_index=index_max_change_alpha;
    end
    if updata_max_change>max_change_u
        updata_max_change=max_change_u;
        updata_index=-index_max_change_u;
    end
    if updata_max_change>max_change_g
        updata_max_change=max_change_g;
        updata_index=index_max_change_g;
    end
    if updata_max_change>max_change_h
        updata_max_change=max_change_h;
        updata_index=-index_max_change_h;
    end
    if updata_max_change>max_change_sum_alpha;
        updata_max_change=max_change_sum_alpha;
        updata_index=0;
    end
end