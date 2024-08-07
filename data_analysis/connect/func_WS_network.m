function [matrix] = func_WS_network(Num,KEtoE,KEtoI,KItoE,KItoI)
% rng();
% clear;clc
% Num=100;
% KEtoE=20;KEtoI=20;KItoE=20;KItoI=20;

PEtoE=KEtoE/100; PEtoI=KEtoI/100; PItoE=KItoE/100; PItoI=KItoI/100;

matrix = zeros(Num);%%%%%
flag_connect=rand(Num,Num);

for ii=1:round(Num)
    for jj=1:round(Num)
        if(ii<=round(Num*0.8)&&jj<=round(Num*0.8)&&flag_connect(ii,jj)<PEtoE)
            matrix(ii,jj)=1;
        end

        if(ii<=round(Num*0.8)&&jj>round(Num*0.8)&&flag_connect(ii,jj)<PEtoI)
            matrix(ii,jj)=1;
        end

        if(ii>round(Num*0.8)&&jj<=round(Num*0.8)&&flag_connect(ii,jj)<PItoE)
            matrix(ii,jj)=1;
        end

        if(ii>round(Num*0.8)&&jj>round(Num*0.8)&&flag_connect(ii,jj)<PItoI)
            matrix(ii,jj)=1;
        end
        if ii==jj, matrix(ii,jj)=0;end
    end
end
% EtoE=sum(matrix(1:round(Num*0.8),1:round(Num*0.8)),'all')/sum(matrix,'all');
% EtoI=sum(matrix(1:round(Num*0.8),round(Num*0.8+1):end),'all')/sum(matrix,'all');
% ItoE=sum(matrix(round(Num*0.8+1):end,1:round(Num*0.8)),'all')/sum(matrix,'all');
% ItoI=sum(matrix(round(Num*0.8+1):end,round(Num*0.8+1):end),'all')/sum(matrix,'all');

EtoE=sum(matrix(1:round(Num*0.8),1:round(Num*0.8)),'all')/80;
EtoI=sum(matrix(1:round(Num*0.8),round(Num*0.8+1):end),'all')/20;
ItoE=sum(matrix(round(Num*0.8+1):end,1:round(Num*0.8)),'all')/80;
ItoI=sum(matrix(round(Num*0.8+1):end,round(Num*0.8+1):end),'all')/20;

% figure (11)
% G = digraph(matrix);
% plot(G)

% G_E=digraph(matrix(1:round(Num*0.8),1:round(Num*0.8)));
% 
% ALL_Cycle=allcycles(G_E);

%% 计算相应的指标
%  [Cc,Cc_avg]          = func_Cluster_Coeff(matrix);
% %  disp(['聚类系数为：',num2str(Cc_avg)]);
%  [Dds_out,Dds_in,Dds_avg_out,Dds_avg_in,P_Dds_out,P_Dds_in]= func_Degree_Distribution(matrix);
% %  disp(['平均度为：',num2str(Dds_avg)]);
%  [Lens,Lens_avg,flag_connect]      = func_Path_Length(matrix);
% %  disp(['平均路径长度为：',num2str(Lens_avg)]);
% 
% 
% 
%  figure (100)
%  subplot(311);
%  bar([1:Num],Dds_out);
%  xlabel('节点编号');
%  ylabel('节点的出度');
% 
%   figure (100)
%  subplot(312);
%  bar([1:Num],Dds_in);
%  xlabel('节点编号');
%  ylabel('节点的入度');

%  subplot(313);
%  bar([0:M],P_Dds,'r');
%  xlabel('节点的度');
%  ylabel('节点度的概率');

% matrix = matrix';