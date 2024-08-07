function [mean_duration,std_duration] = fun_std(T_persistent)

% clear;clc
% load('matrxi_15x100_0.03.mat')
N=size(T_persistent,1);
std_duration(1:N)=0;

mean_duration = mean (T_persistent,2);
for ii=1:N
    std_duration(ii) = std (T_persistent(ii,:));
end
