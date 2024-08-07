function [Dds_out,Dds_out_to_E,Dds_out_to_I,Dds_in,Dds_in_from_E,Dds_in_from_I,Dds_avg_out,Dds_avg_in,P_Dds_out,P_Dds_in]=func_Degree_Distribution(matrix)
 
Num = size(matrix,2);

%% 出度
Dds_out = zeros(1,Num);Dds_out_to_E= zeros(1,Num);Dds_out_to_I= zeros(1,Num);
for i=1:Num
    Dds_out(i)=sum(matrix(i,:));
    Dds_out_to_E(i)=sum(matrix(i,1:round(0.8*Num)));
    Dds_out_to_I(i)=sum(matrix(i,round(0.8*Num+1):Num));
end
Dds_avg_out = mean(Dds_out);

M=max(Dds_out);
for i=1:M+1
    Num_Dds_out(i) = length(find(Dds_out==i-1));
end
P_Dds_out    = zeros(1,M+1);
P_Dds_out(:) = Num_Dds_out(:)./sum(Num_Dds_out);
%% 入度
Dds_in = zeros(1,Num);Dds_in_from_E= zeros(1,Num);Dds_out_to_I= zeros(1,Num);
for i=1:Num
    Dds_in(i)=sum(matrix(:,i));
    Dds_in_from_E(i)=sum(matrix(1:round(0.8*Num),i));
    Dds_in_from_I(i)=sum(matrix(round(0.8*Num+1):Num,i));
end
Dds_avg_in = mean(Dds_in);

M=max(Dds_in);
for i=1:M+1
    Num_Dds_in(i) = length(find(Dds_in==i-1));
end
P_Dds_in    = zeros(1,M+1);
P_Dds_in(:) = Num_Dds_in(:)./sum(Num_Dds_in);

 



