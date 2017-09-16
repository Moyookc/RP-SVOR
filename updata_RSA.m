function [out_R]=updata_RSA(x,svm_struct,index,R,y)
support_set=svm_struct.index_S_S;
if index>0  
    local_index1=find(support_set==index);
    if length(local_index1)==0                         %increment
        R=increment_RSA(R,index,x,svm_struct,y);
    else                                               %decrement
        pre_length=1+local_index1;
        R=decrement_R(R,pre_length);
    end
end
out_R=R
tmp=initialR(x,y,svm_struct);