function [x,y,Q,c,q,Aeq,beq,peq,Aiq,biq,piq,L,U,O]=PrepareMatsORwei(t_range,flag)
global KTYPE KSCALE  level size_training
%     KTYPE =6;
%     KSCALE =0.5;
    start_u=0.1;
    epsilon=1;
    [x,y]=read_data(flag);
    level=max(y);
    x=zscore(x);
    [x,y,test_x,test_y]=random_select(x,y,size_training);
    [local_x,local_y,length_data]=trans_data(x,y);
    local_length=ones(1,level-1)*length_data(1:level-1)+ones(1,level-1)*length_data(2:level);
    Q=kernel(local_x,local_x).*(local_y*local_y');
    tmp1=zeros(local_length,level-2);
    local_length2=local_length+level-2;
    Q=[Q  tmp1];
    tmp2=zeros(level-2,local_length2);
    Q=[Q;tmp2];
    c=[-ones(local_length,1);zeros(level-2,1)];
    q=zeros(local_length2,1);
    [Aeq,beq,peq,Aiq,biq,piq]=getAeqwei(length_data,local_length);
    L=zeros(local_length2,1);
    U=[ones(local_length,1)*start_u;ones(level-2,1)*inf];
    O=[ones(local_length,1);zeros(level-2,1)];
end
function [Aeq,beq,peq,Aiq,biq,piq]=getAeqwei(length_data,sum_length)
    global level;
    local_length2=sum_length+level-2;
    Aeq=zeros((level-1),local_length2);
    Aiq=zeros((level-1),local_length2);
    index1=0;
    index2=ones(1,level-1)*length_data(1:level-1);
    for i=1:(level-1)
        Aeq(i,index1+1:index1+length_data(i))=-ones(1,length_data(i));
        Aeq(i,index2+1:index2+length_data(i+1))=ones(1,length_data(i+1));
        index1=index1+length_data(i);
        index2=index2+length_data(i+1);
        if i==1
            Aeq(i,sum_length+i)=1;
        else
            if i==(level-1)
                Aeq(i,sum_length+i-1)=-1;
            else
                Aeq(i,sum_length+i)=1;
                Aeq(i,sum_length+i-1)=-1;
            end
        end
    end
    Aeq=Aeq';
    Aiq=[];
    beq=zeros((level-1),1);
    peq=zeros((level-1),1);
    biq=[];
    piq=[];
end