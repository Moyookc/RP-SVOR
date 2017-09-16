function out_KKT=sub_KKT2(alpha,local_g,svm_struct)
    global fake_zero size_training
    out_KKT=1;
    C=1/size_training;
    index_S_S=svm_struct.index_S_S;
    sum_length=length(index_S_S);
    for i=1:sum_length
        i1=index_S_S(i);
        if local_g(i1)>fake_zero | local_g(i1)<-fake_zero
            out_KKT=0;
            break;
        end                    
    end
end