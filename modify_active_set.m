function [svm_struct,g,Qc]=modify_active_set(svm_struct,g,Qc,y)
    index_S_S=svm_struct.index_S_S;
    y_S_S=y(index_S_S);
    svm_struct=svm_struct;
    g=g;
    Qc=Qc;
    if sum(y_S_S==1)==0
        
    elseif sum(y_S_S==-1)==0

    end    
end