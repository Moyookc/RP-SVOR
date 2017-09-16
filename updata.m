function    [out_g,out_h,out_struct,out_Qc]=updata...
        (mystruct,g,h,beta,gama_g, gama_h,max_delta_zeta,index,Q,q,Aeq,peq,Aiq,piq,L,U,O)
    [out_struct]=updata_struct(mystruct,beta,max_delta_zeta,index,L,U,O);
    [out_g,out_h]=updata_gh(g,h,mystruct,gama_g, gama_h,max_delta_zeta);
%    [out_svm_struct,out_g]=Adjust4EmptySS(out_svm_struct,out_g,x,y);
    out_Qc=getQcSA(out_struct,Q,q,Aeq,peq,Aiq,piq,O);
end