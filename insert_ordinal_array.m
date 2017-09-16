function [out,ordinal]=insert_ordinal_array(in,index)
    local_length=length(in);
    ordinal=1:local_length;
    ordinal=ordinal';
    local_tmp=find(in<index);
    local_length1=length(local_tmp);
    if local_length1<local_length & local_length1>0
        in=[in(1:local_length1);index;in(local_length1+1:local_length)];
        ordinal=[ordinal(1:local_length1); local_length+1; ordinal(local_length1+1:local_length)];
    else
        if local_length1==0
            in=[index;in];
            ordinal=[local_length+1;ordinal];
        end
        if local_length1==local_length
            in=[in;index];
            ordinal=[ordinal;local_length+1];
        end
    end
    out=in;
end