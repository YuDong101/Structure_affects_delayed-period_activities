function [Dds_EfromE,Dds_EfromI,Dds_IfromE,Dds_IfromI]=func_Degree_Distribution_E_from_I(matrix)
 
Num = size(matrix,2);
%% 出度
Dds_EfromE = zeros(1,round(Num*0.8));
for i=1:round(Num*0.8)
    Dds_EfromE(i)=sum(matrix(1:round(Num*0.8),i));
end

Dds_EfromI = zeros(1,round(Num*0.2));
for i=1:round(Num*0.8)
    Dds_EfromI(i)=sum(matrix(round(Num*0.8+1):Num,i));
end

Dds_IfromI = zeros(1,round(Num*0.2));
for i=round(Num*0.8+1):Num
    Dds_IfromI(i)=sum(matrix(round(Num*0.8+1):Num,i));
end

Dds_IfromE = zeros(1,round(Num*0.2));
for i=round(Num*0.8+1):Num
    Dds_IfromE(i)=sum(matrix(1:round(Num*0.8),i));
end

% Dds_avg_out = mean(Dds_fromE);
% 
% M=max(Dds_fromE);
% for i=1:M+1
%     Num_Dds_out(i) = length(find(Dds_fromE==i-1));
% end
% P_Dds_out    = zeros(1,M+1);
% P_Dds_out(:) = Num_Dds_out(:)./sum(Num_Dds_out);