function  [Cc,Cc_avg]=func_Cluster_Coeff(matrix)

Num = size(matrix,2);
Cc  = zeros(1,Num);

for i=1:Num
    Index1 = find(matrix(i,:)==1);  
    if isempty(Index1) == 1
       Cc(i)=0;
    else
        m = length(Index1); 
        if m == 1
           Cc(i)=0;
        else
           B     = matrix(Index1,Index1);   
           Cc(i) = length(find(B==1))/(m*(m-1));
        end
    end
end
Cc_avg=mean(Cc);

% [matrix] = matrix_trans(matrix_rec_R,1,1);
% 
% N = size(matrix,2);
% D1=zeros(N,N);
% for i=1:N
%     for j=1:N
%     D1(i,j)=matrix(i,j);
%     end
% end
% %超边度的矩阵
% node_degree=zeros(1,N);
% for h=1:N
%  node_degree(1,h)=sum(D1(h,:));
% end
%  node_degree;
%  cc=zeros(1,N);
%  %求聚类系数
%  for i=1:N
%  l1=0;
%  s=zeros(1,N);
%  if node_degree(1,i)>=2
%      for j=1:N
%           if D1(i,j)==1
%              s(1,j)=1;
%           end
%      end
%      
%      for k=1:N
%        if s(1,k)==1
%            s1=k;
%            for x=s1:N
%                if s(1,x)==1
%                     s2=x;
%                     if D1(s1,s2)==1
%                       l1=l1+1;
%                     end
%                end
%            end
%        end
%      end
%      cc(1,i)=(2*l1)/((node_degree(1,i)*(node_degree(1,i)+1)));
%  end
%  end
%  
%  sum1=0;
%  for i=1:N
%   sum1=sum1+cc(1,i)*1;
%  end
% 
% %平均聚类系数
%  ave=sum1/N;
% 
% 

% A=matrix;
% N=size(A,2);
% D=0;
% for i=1:N
%     C=[];
%     for j=1:N
%         if A(i,j)==1
%             C=[C,j];
%         end
%     end
%     K=size(C,2);
%     B=zeros(K);
%     for p=1:K
%         for q=p:K
%             B(p,q)=A(C(1,p),C(1,q));
%         end
%     end
%     D=D+(sum(B(:)))/(K*(K-1)/2);
% end
% y=D/N;

