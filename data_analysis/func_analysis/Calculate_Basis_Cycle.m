function [sum_Num_Cycle,sum_Num_Cycle_E] = Calculate_Basis_Cycle(matrix)

% clear;clc
% load matrix.mat

[Num,uesless]=size(matrix);
G = digraph(matrix);
G_E=digraph(matrix(1:round(Num*0.8),1:round(Num*0.8)));

% ALL_Cycle=allcycles(G,'MinCycleLength',3);  toc % ,'MinCycleLength',3,'MaxCycleLength',6

matrix_ul=triu(matrix)+tril(matrix)';
matrix_ul(matrix_ul~=0)=1;
G_C = graph(matrix_ul,"upper");

cycle_basis=cyclebasis(G_C);   % 计算无向图的基循环

[Num_cyclebasis,useless]=size(cycle_basis);
cycle_basis_dir={}; cycle_basis_dir_E={};
for ii=1:Num_cyclebasis   % 由邻接矩阵限制得到有向图的基循环
    arrayCycle=cycle_basis{ii};
    lll=length(arrayCycle);
    for jj=1:lll+1
        if jj==lll+1, 
            aaa=find(arrayCycle>0.8*Num);
            if isempty(aaa),cycle_basis_dir_E{ii,1}=arrayCycle;
            else cycle_basis_dir{ii,1}=arrayCycle;end
            break;
        end
        if jj<lll && matrix(arrayCycle(jj),arrayCycle(jj+1))==0, break; end
        if jj==lll && matrix(arrayCycle(jj),arrayCycle(1))==0, break; end
    end
end

cycle_basis_dir(cellfun(@isempty,cycle_basis_dir))=[];  % 清除为空的元素
cycle_basis_dir_E(cellfun(@isempty,cycle_basis_dir_E))=[];  % 清除为空的元素
[Num_cdir,useless]=size(cycle_basis_dir);
[Num_cdir_E,useless]=size(cycle_basis_dir_E);

Num_cyclebasis_dir=[]; Num_cyclebasis_dir_E=[];
for iii=1:Num_cdir    % 得到基循环的边长
    Num_cyclebasis_dir(iii)=length(cycle_basis_dir{iii});
end
for iii=1:Num_cdir_E    % 得到基循环的边长
    Num_cyclebasis_dir_E(iii)=length(cycle_basis_dir_E{iii});
end

Num_Cycle_max=max(Num_cyclebasis_dir); Num_Cycle_min=min(Num_cyclebasis_dir);
Num_Cycle_max_E=max(Num_cyclebasis_dir_E); Num_Cycle_min_E=min(Num_cyclebasis_dir_E);

sum_Num_Cycle(1:10)=0.0;
for jjj=Num_Cycle_min:Num_Cycle_max
    sum_Num_Cycle(jjj)=sum(Num_cyclebasis_dir==jjj);
end

sum_Num_Cycle_E(1:10)=0.0;
for jjj=Num_Cycle_min_E:Num_Cycle_max_E
    sum_Num_Cycle_E(jjj)=sum(Num_cyclebasis_dir_E==jjj);
end

% P_Num_Cycle=sum_Num_Cycle/sum(sum_Num_Cycle);



