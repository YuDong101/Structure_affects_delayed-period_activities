function [x_E,output_E,x_I,output_I,Cc_avg,Lens_avg] = matrix_analysis(matrix,K)
% clear;clc; rng();
% % load mar.mat
% Num=5; K=2; Per=1.0;
% [matrix] = func_WS_network(Num,K,Per);

M = matrix';
n = size(M,1);  x_E(1,n+1)=0; output_E(1,n+1)=0; x_I(1,n+1)=0; output_I(1,n+1)=0;
W = M; N_loop=0;
%% 计算相应的指标
 [Cc,Cc_avg]          = func_Cluster_Coeff(matrix);
%  disp(['聚类系数为：',num2str(Cc_avg)]);
 [Dds_out,Dds_in,Dds_avg_out,Dds_avg_in,P_Dds_out,P_Dds_in]= func_Degree_Distribution(matrix);
%  disp(['平均度为：',num2str(Dds_avg)]);
 [Lens,Lens_avg,flag_connect]      = func_Path_Length(matrix);
%  disp(['平均路径长度为：',num2str(Lens_avg)]);
if flag_connect==0,return; end
 %% 分成子环
for jjj=1:n
    k = zeros(1,n+1); a = zeros(n+1,K);
    v = jjj; k(1) = v;
    s=1; i=2; j=1; flag_runLoop=1;
    while i<=n+1
%          s
    %         fprintf('flag_runLoop=%d\n',flag_runLoop);
        if flag_runLoop==1, a(i,:) = find(M(k(i-1),:) == 1);end
        flag_runLoop=0;
        non_zero = find(a(i,:)~=0,1); % 找到一行中第一个非0的元素
        k(i:end)=0;
%         a
    %         fprintf('i=%d non_zero=%d \n',i,non_zero);
        
        if isempty(non_zero) && K == find(a(i-1,:)~=0,1)
            i=i-1;
            non_zero = find(a(i-1,:)~=0,1);
            a(i-1,non_zero) = 0;
            non_zero = find(a(i-1,:)~=0,1);
            k(i-1)=a(i-1,non_zero);
            if isempty(find(k(2:end)==a(i-1,non_zero),K)), flag_runLoop=1;
            else, i=i-1;
            end
        end
%             k
        
        if isempty(non_zero) &&  K > find(a(i-1,:)~=0,1)
    %             i=i-1;
            non_zero = find(a(i-1,:)~=0,1); % 这四行标配
            a(i-1,non_zero) = 0;
            non_zero = find(a(i-1,:)~=0,1);
            k(i-1)=a(i-1,non_zero);
            if isempty(find(k(2:end)==a(i-1,non_zero),K)), flag_runLoop=1;
            else, i=i-1;
            end
        end
%         if s==415, break;end
    
        if flag_runLoop==1, a(i,:) = find(M(k(i-1),:) == 1);end
        flag_runLoop=1;
        non_zero = find(a(i,:)~=0,1); % 找到一行中第一个非0的元素
        k(i:end)=0;
        
        if size(non_zero,2)>=1 && isempty(find(k(2:end)==a(i,non_zero),K))  % a(i,1) 不属于 k(2:end)
            k(i)=a(i,non_zero);    % 记录环路顺序
        else%if size(non_zero,2)>=1
            k(i)=a(i,non_zero);    % 记录环路顺序
            a(i,non_zero) = 0;%  'Node recurrence'
            i=i-1;
            flag_runLoop=0;
        end

        non_zero = find(a(i,:)~=0,1); % 找到一行中第一个非0的元素
        if a(i,non_zero)==k(1) % 第一个非0元素已经返回 v(顶点)，回路结束。回退上一步，找其他回路
            a(i,non_zero) = 0;%  'Loop record'
            i = i-1;
            j=j+1;
            k_loop(jjj,j,:)=k;
            N_loop(jjj,j) = i;          % 只有返回 v 才记录
            flag_runLoop=0;
        end
        a_revised=a;
        
    %         k
    %         i
%         ii=i
        if i>2.1 && i<=n, non_zero_revised = find(a(i+1,:)~=0,1);non_zero = find(a(i,:)~=0,1);end
        if i>2.1 && i<=n && flag_runLoop==0 && isempty(non_zero_revised) && K==non_zero
%             ss=0;
            while  K == find(a(i-1,:)~=0,1)
                if K == find(a(i,:)~=0,1) && K == find(a(i-1,:)~=0,1)                    
                    a(i,K) = 0; 
                    i=i-1;
%                      'yes'
                end
%                 'yyes'
%                  if ss==20, break;end
%                  ss=ss+1;
            end
        end

    %         a_revised=a
%         if s==420, break;end
        i=i+1;
        s=s+1;
        if isempty(find(a(3,:)~=0,1))==1 && find(a(2,:)~=0,1)==K && s>2, break; end
    end
end

n_k=j;
for jjj=1:n
    for j=1:n_k
        for iii=round(0.8*n+1):n
            if isempty(find(k_loop(jjj,j,:)==iii,1))  %非空
                N_loop_E(jjj,j) = N_loop(jjj,j);
            else
                N_loop_I(jjj,j) = N_loop(jjj,j);
            end
        end
    end
end
 [x_out,output_out] = fun_histogram(N_loop);
 [x_E,output_E] = fun_histogram(N_loop_E);
 [x_I,output_I] = fun_histogram(N_loop_I);

%  figure (100)
%  subplot(211);
%  bar([1:n],Dds_out);
%  xlabel('节点编号');
%  ylabel('节点的出度');
% 
%   figure (100)
%  subplot(212);
%  bar([1:n],Dds_in);
%  xlabel('节点编号');
%  ylabel('节点的入度');

%  subplot(212);
%  bar([0:M],P_Dds,'r');
%  xlabel('节点的度');
%  ylabel('节点度的概率');



% subplot (2,1,1), plot (x_E,output_E);
% subplot (2,1,2), plot (x_I,output_I);
    