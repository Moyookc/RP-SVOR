function out_KKT=test_KKT2(svm_struct,g,y)
global fake_zero 
out_KKT=0;
local_g=g;
alpha=svm_struct.alpha;
local_length=length(y);
zero=alpha'*y;
if zero<fake_zero & zero>-fake_zero
    sum=alpha'*ones(local_length,1);
    if sum>svm_struct.nu-fake_zero & sum<svm_struct.nu+fake_zero     
        local_tmp=sub_KKT(alpha,local_g);
        if local_tmp==1
            out_KKT=1;
        end
    end
end