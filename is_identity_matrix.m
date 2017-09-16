function flag=is_identity_matrix(R)
global fake_zero
local_length=length(R);
flag=1;
for i=1:local_length
    for j=1:local_length
        a=R(i,j);
        if i==j
            if a>1+fake_zero | a<1-fake_zero
                flag=0;
                break;
            end
        else
            if a>fake_zero | a<-fake_zero
                flag=0;
                break;
            end
        end        
    end
end