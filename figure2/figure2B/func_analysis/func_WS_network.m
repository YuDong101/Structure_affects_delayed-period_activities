function [matrix] = func_WS_network(Num,K,Per)
% rng(1);
%����С����
angle = 2*pi/Num:2*pi/Num:2*pi;  
x     = 100*sin(angle);
y     = 100*cos(angle);
matrix= zeros(Num);
%������������
for i1=1:Num
    for j1=i1+1:i1+K
        j2=j1;
        if j1 > Num
           j2 = mod(j1,Num);
        end
      matrix(i1,j2) = 1;
%       matrix(j2,i1) = 1;  ????
    end
end

for i1=1:Num
    for j1=i1+1:i1+K
        j2=j1;
        if j1>Num
           j2=mod(j1,Num);
        end
        p1 = rand();
        %������������ж��Ƿ��������
        if p1 < Per             
           matrix(i1,j2) = 0;
%            matrix(j2,i1) = 0;
           matrix(i1,i1) = inf;
           a             = find(matrix(i1,:)==0);
           rand_data     = randi([1,length(a)],1,1);
           jjj           = a(rand_data);
           matrix(i1,jjj)= 1;
%            matrix(jjj,i1)= 1;
           matrix(i1,i1) = 0;           
        end
    end
end

matrix = matrix';