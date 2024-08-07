function [x,y]=trans_fitResults(fitResult,x_data)
% fitResult=fitResult_All_Cycle_t;
% x_data=x_All_Cycle;

N_point=101;
x=min(x_data):(max(x_data)-min(x_data))/(N_point-1):max(x_data);

for ii=1:N_point
    y(ii) = fitResult.a1 * exp(-((x(ii) - fitResult.b1)/fitResult.c1)^2);
end