clc;
clear;
close all;
warning off;
addpath 'func\'
rng(1);

%�����������
%Num--���������Per--���Ӹ��ʣ�avg_length--�ߵ�ƽ������
Num          = 100;
K            = 3;
Per          = 0.4;
[matrix,x,y] = func_WS_network(Num,K,Per);

%������efficiency
Pset         = [0.99:-0.001:0.01];
Eg           = [];
for kk = 1:length(Pset)
    kk
    matrixs = matrix;
    for i=1:Num 
        for j=i+1:Num
            if matrix(i,j)~=0
               %��һ�����ʽ��й���
               p = rand;
               if p < Pset(kk)
                  matrixs(i,j) = 0; 
               end
            end
        end
    end
    %����efficiency
    Es = 0;
    for i=1:Num 
        for j=i+1:Num
            if matrixs(i,j)~=0 & (sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2)) > 2000*Pset(kk)
               Es = Es + 3e3/(sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2)); 
            end
        end
    end
    Eg = [Eg,Es/Num/(Num-1)];
end 

figure;
plot(Pset,Eg,'linewidth',2);
xlabel('p = fraction of removed nodes'); 
ylabel('efficiency'); 
save atrack3.mat Pset Eg
axis([0,1,0,0.5]);



