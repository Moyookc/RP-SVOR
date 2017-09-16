function out_R=decrement_R(R,pre_length)    

    local_tmp1=R(pre_length,:);
    local_tmp2=R(:,pre_length);
    local_xx=R(pre_length,pre_length);
    local=local_tmp2*local_tmp1/local_xx;
    out_R=R-local;
    out_R(pre_length,:)=[];
    out_R(:,pre_length)=[];