function [sum_Num_Cycle,sum_Num_Cycle_E] = Calculate_All_Cycle(matrix)

% clear;clc
% load matrix.mat

[Num,uesless]=size(matrix);
G = digraph(matrix);
G_E=digraph(matrix(1:round(Num*0.8),1:round(Num*0.8)));

ALL_Cycle=allcycles(G,'MinCycleLength',3,'MaxCycleLength',3); % ,'MinCycleLength',3,'MaxCycleLength',6

[Num_cyclebasis,useless]=size(ALL_Cycle);
ALL_Cycle_dir={}; ALL_Cycle_dir_E={};
for ii=1:Num_cyclebasis   % 由邻接矩阵限制得到有向图的基循环
    arrayCycle=ALL_Cycle{ii};
    lll=length(arrayCycle);
    for jj=1:lll+1
        if jj==lll+1, 
            aaa=find(arrayCycle>0.8*Num);
            if isempty(aaa),ALL_Cycle_dir_E{ii,1}=arrayCycle;
            else ALL_Cycle_dir{ii,1}=arrayCycle;end
            break;
        end
        if jj<lll && matrix(arrayCycle(jj),arrayCycle(jj+1))==0, break; end
        if jj==lll && matrix(arrayCycle(jj),arrayCycle(1))==0, break; end
    end
end

ALL_Cycle_dir(cellfun(@isempty,ALL_Cycle_dir))=[];  % 清除为空的元素
ALL_Cycle_dir_E(cellfun(@isempty,ALL_Cycle_dir_E))=[];  % 清除为空的元素
[Num_cdir,useless]=size(ALL_Cycle_dir);
[Num_cdir_E,useless]=size(ALL_Cycle_dir_E);

Num_ALL_Cycle_dir=[]; Num_ALL_Cycle_dir_E=[];
for iii=1:Num_cdir    % 得到基循环的边长
    Num_ALL_Cycle_dir(iii)=length(ALL_Cycle_dir{iii});
end
for iii=1:Num_cdir_E    % 得到基循环的边长
    Num_ALL_Cycle_dir_E(iii)=length(ALL_Cycle_dir_E{iii});
end

Num_Cycle_max=max(Num_ALL_Cycle_dir); Num_Cycle_min=min(Num_ALL_Cycle_dir);
Num_Cycle_max_E=max(Num_ALL_Cycle_dir_E); Num_Cycle_min_E=min(Num_ALL_Cycle_dir_E);

sum_Num_Cycle(1:10)=0.0;
for jjj=Num_Cycle_min:Num_Cycle_max
    sum_Num_Cycle(jjj)=sum(Num_ALL_Cycle_dir==jjj);
end

sum_Num_Cycle_E(1:10)=0.0;
for jjj=Num_Cycle_min_E:Num_Cycle_max_E
    sum_Num_Cycle_E(jjj)=sum(Num_ALL_Cycle_dir_E==jjj);
end

% P_Num_Cycle=sum_Num_Cycle/sum(sum_Num_Cycle);



