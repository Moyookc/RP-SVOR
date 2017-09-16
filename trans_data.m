function [local_x,local_y,length_data]=trans_data(training_x,training_y);
    global level;
    local_x=[];
    for i=1:level-1
        local_tmp=find(training_y==i);
        local_x=[local_x;training_x(local_tmp,:)];
    end
    for i=2:level
        local_tmp=find(training_y==i);
        local_x=[local_x;training_x(local_tmp,:)];
    end

    length_data=zeros(level,1);
    for i=1:level;
        length_data(i)=length(find(training_y==i));
    end
    local_y=[-ones(sum(length_data([1:level-1])),1);ones(sum(length_data([2:level])),1)];
end