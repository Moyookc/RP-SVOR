function out_flag=test_R(x,y,svm_struct,R)
global epsilon
index_S_S=svm_struct.index_S_S;
length_S_S=length(index_S_S);
x_S_S=x(index_S_S,:);
y_S_S=y(index_S_S);
Q=kernel(x_S_S,x_S_S).*(y_S_S*y_S_S');
Q=[y_S_S';Q];
Q=[Q;ones(1,length_S_S)];
lable1=[0;ones(length_S_S,1);epsilon];
label2=[0;y_S_S;0];
Q=[label2 Q];
Q=[ Q lable1];
out=Q*R;
out_flag=is_identity_matrix(out);
