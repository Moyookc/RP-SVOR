function [x,y]=read_data(flag)
    if flag==1
        load Bank;
    end
    if flag==2
        load Computer;
    end
    if flag==3
        load Abalone;
    end
    if flag==4
        load Wine_red;
    end
    if flag==5
        load Wine_white;
    end
    if flag==6
        load Spine;
    end
    x=zscore(x);
end
