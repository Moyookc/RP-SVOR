function flag=ChekOneClass(svm_struct,y)
support_set=svm_struct.index_S_S;
local_length=length(support_set);
local_y=y(support_set);
sum=local_y'*ones(local_length,1);
sum=abs(sum);
flag=0;
if sum==local_length
    flag=1;
end