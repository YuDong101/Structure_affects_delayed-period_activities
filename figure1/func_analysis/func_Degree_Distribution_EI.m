function [Dds_fromE,Dds_fromI]=func_Degree_Distribution_EI(matrix)
 
Num = size(matrix,2);
%% 出度
Dds_fromE = zeros(1,round(Num*0.8));
for i=1:Num
    Dds_fromE(i)=sum(matrix(1:round(Num*0.8),i));
end

Dds_fromI = zeros(1,round(Num*0.2));
for i=1:Num
    Dds_fromI(i)=sum(matrix(round(Num*0.8+1):Num,i));
end

% Dds_avg_out = mean(Dds_fromE);
% 
% M=max(Dds_fromE);
% for i=1:M+1
%     Num_Dds_out(i) = length(find(Dds_fromE==i-1));
% end
% P_Dds_out    = zeros(1,M+1);
% P_Dds_out(:) = Num_Dds_out(:)./sum(Num_Dds_out);