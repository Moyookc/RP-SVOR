function [out_struct]=updata_svm_struct(mystruct,...
    beta,max_delta_t,index,L,U,O)
    global epsilon fake_zero 
    alpha=mystruct.alpha;
    SE=mystruct.SE;
    SR=mystruct.SR;
    len1=length(SE);
    beta_alpha=beta(1:len1);
    alpha(SE)=alpha(SE)+beta_alpha*max_delta_t;
    alpha(SR)=alpha(SR)+O(SR)*max_delta_t;
    mystruct.alpha=alpha;
    lang_u=mystruct.lang_u;
    indexJ=mystruct.indexJ;
    len2=length(indexJ);
    if len2>0
        beta_u=beta(len1+1:end);
        lang_u(indexJ)=lang_u(indexJ)+beta_u*max_delta_t;
        mystruct.lang_u=lang_u;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%t
    t=mystruct.t;
    t=t+max_delta_t;
    mystruct.t=t;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%set
    SE=mystruct.SE;
    SL=mystruct.SL;
    SR=mystruct.SR;
    indexCJ=mystruct.indexCJ;
    if index>0
        flag_index=inf;
        if alpha(index)>U(index)+O(index)*t-fake_zero % SE
            flag_index=0;
        end
        if alpha(index)<L(index)+fake_zero % SR
            flag_index=1;
        end
        local_index_SE=find(SE==index);
        local_index_SL=find(SL==index);
        local_index_SR=find(SR==index);
        if length(local_index_SE)>0
            SE(local_index_SE)=[];
            if flag_index==0
                [SR,ordinal]=insert_ordinal_array(SR,index);
            end
            if flag_index==1
                [SL,ordinal]=insert_ordinal_array(SL,index);
            end
        else
            if length(local_index_SL)>0
                SL(local_index_SL)=[];
                [SE,ordinal]=insert_ordinal_array(SE,index);
            else
                SR(local_index_SR)=[];
                [SE,ordinal]=insert_ordinal_array(SE,index);
            end
        end
    end
    if index<0
        local_index_J=find(indexJ==index);
        local_index_CJ=find(indexCJ==index);
        if length(local_index_J)>0
            indexJ(local_index_J)=[];
            [indexCJ,ordinal]=insert_ordinal_array(indexCJ,index);
        else
            indexCJ(local_index_CJ)=[];
            [indexJ,ordinal]=insert_ordinal_array(indexJ,index);
        end
    end
    mystruct.SE=SE;
    mystruct.SR=SR;
    mystruct.SL=SL;
    mystruct.indexJ=indexJ;
    mystruct.indexCJ=indexCJ;
    out_struct=mystruct;
end
