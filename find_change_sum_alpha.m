function [max_change_sum_alpha]=find_change_sum_alpha(mystruct,t_range)
    t=mystruct.t;
    max_change_sum_alpha=t_range(2)-t;
end