
function [matrix] = func_RtoS(NumR,NumS,Per)
% clear;clc
% NumS=10; NumR=10; Per=0.3;
% rng(1);
K=round(Per*NumS);
matrix = zeros(round(NumS*0.8),NumR);%%%%%

for j1=1:NumR
    for ii=1:K
        matrix(j1,j1) = inf;
        a             = find(matrix(j1,1:round(NumS*0.8))==0);
        rand_data     = randi([1,length(a)],1,1);
        jjj           = a(rand_data);
        matrix(j1,jjj)= 1;
        matrix(j1,j1) = 0;
    end
end

% matrix(round(NumR*0.8+1):NumR,NumS) = 0;

matrix = matrix';