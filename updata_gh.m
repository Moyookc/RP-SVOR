function [out_g,out_h]=updata_gh(g,h,mystruct,gama_g, gama_h,max_delta_zeta)
    other_set=[mystruct.SR;mystruct.SL];
    g(other_set)=g(other_set)+gama_g(other_set)*max_delta_zeta;
    out_g=g;
    indexCJ=mystruct.indexCJ;
    h(indexCJ)=h(indexCJ)+gama_h(indexCJ)*max_delta_zeta;
    out_h=h;
end