function FlagUnworkable=CheckUnworkable(svm_struct,g,y)
global fake_zero
FlagUnworkable=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%find S_E
other_set=svm_struct.index_other;
local_length=length(svm_struct.alpha);
label_other=zeros(local_length,1);
label_other(other_set)=1;
g_noless_zero=g<-fake_zero;
label_S_R=label_other.*g_noless_zero;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%find S_{-y}
support_set=svm_struct.index_S_S;
y_S_S=y(support_set(1));
label_negative_y=(y==-y_S_S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
local_tmp=label_S_R.*label_negative_y;
local_tmp2=length(find(local_tmp==1));
if local_tmp2==0
    FlagUnworkable=1;
end


