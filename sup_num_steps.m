clear all;
global KTYPE size_training KSCALE solver_type
solver_type=2;   % 1: quadprog, 2: smo solver
data_flag=4;
kscale=[1 0.01 0.01 0.1 1];
KTYPE=6
KSCALE=0.1
size_training=300;
[out]=main(data_flag);

