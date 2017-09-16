
% RP-SVOR (Regularization Path Algorithm for Support Vector Ordinal Regression) demo 
% 
% Bin Gu, Nanjing Univerity of Information Secience and Technolgy
% 2016.3.23
% 
clear all;
global KTYPE size_training KSCALE solver_type
solver_type=2;   % 1: quadprog, 2: smo solver
data_flag=4;  % 1-6
KTYPE=6  % kernel type: 1: linear, 2, polynomal, 6: gaussian
KSCALE=1 % paramter in kernel, please find details in kernel.m
size_training=400; % size of dataset by randomly selecting
[out]=main(data_flag);  % [out] is the solution path returned by our RP-SVOR
