
function [matrix] = func_StoR(NumS,NumR,Per)
% clear;clc
% NumS=10; NumR=10; Per=0.3;
% rng(2);
K=round(Per*NumS);
matrix = zeros(NumS,round(0.8*NumR));%%%%%

for j1=1:round(NumR*0.8)
    for ii=1:K
        matrix(j1,j1) = inf;
        a             = find(matrix(j1,:)==0);
        rand_data     = randi([1,length(a)],1,1);
        jjj           = a(rand_data);
        matrix(j1,jjj)= 1;
        matrix(j1,j1) = 0;
    end
end

matrix(NumS,round(NumR*0.8+1):NumR) = 0;

matrix = matrix';



