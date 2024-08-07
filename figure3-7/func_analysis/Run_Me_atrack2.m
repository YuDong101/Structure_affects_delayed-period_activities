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

%�������,ģ��500�Σ������ܿ�������
X1 = [];
X2 = [];
X3 = [];

for kk = 1:500
    kk
    matrixs               = matrix;
    [Dds,Dds_avg,M,P_Dds] = func_Degree_Distribution(matrixs);  
    %���ݽڵ�ȴ�С����ĳһ��Χ�ڵĶȵĽڵ���й���
    Dmax = max(Dds);
    ind2 = find(Dds>=0.7*Dmax);
    for i1 = 1:length(ind2)
        for i=1:Num 
            if i == ind2(i1)
                for j=i+1:Num
                    if matrix(i,j)~=0
                       %��һ�����ʽ��й���
                       p = rand;
                       if p > 0.9
                          matrixs(i,j) = 0; 
                       end
                    end
                end
            end
        end
    end

    %������Ӧ��ָ��
    [Cc,Cc_avg]          = func_Cluster_Coeff(matrixs);
    [Dds,Dds_avg,M,P_Dds]= func_Degree_Distribution(matrixs);  
    [Lens,Lens_avg]      = func_Path_Length(matrixs);   
    if Lens_avg < 10000
       X1 = [X1,Cc_avg];
       X2 = [X2,Dds_avg];
       X3 = [X3,Lens_avg];
    end
end

figure;
subplot(311);
plot(X1);
xlabel('ģ�������������');
ylabel('����ϵ���仯');
subplot(312);
plot(X2);
xlabel('ģ�������������');
ylabel('ƽ���ȱ仯'); 
subplot(313);
plot(X3);
xlabel('ģ�������������');
ylabel('ƽ��·�����ȱ仯');  
 
save atrack2.mat X1 X2 X3 
 
std(X1/max(X1))
std(X2/max(X2))
std(X3/max(X3))








