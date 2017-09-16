function [training_x,training_y,test_x,test_y]=random_select(x,y,l);
global level;
local_flag=1;
while(local_flag)
    training_x=[];
    training_y=[];
    test_x=[];
    test_y=[];
%    local_index=L;
    local_length=length(x(:,1));
    test_x=x;
    test_y=y;
    for i=1:l
        local_rand= randint(1,1,[1,local_length+1-i]);
        local_x=test_x(local_rand,:);
        local_y=test_y(local_rand);
        training_x=[training_x;local_x];
        training_y=[training_y;local_y];
        test_x(local_rand,:)=[];
        test_y(local_rand)=[];        
    end
    for i=1:level
        local_flag=local_flag*length(find(training_y==i));
    end
    if local_flag==0
        local_flag=1;
    else
        local_flag=0;
    end
end