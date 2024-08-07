function [x_tran,y_tran,fitResult] = gaussfitt(x,y)
% x=xCc;y=hp_Cc;
% 定义高斯模型
gaussianModel = fittype('a1 * exp(-((x - b1)/c1)^2)', 'independent', 'x', 'dependent', 'y');
% 进行拟合
fitResult = fit(x', y', gaussianModel, 'StartPoint', [max(y), mean(x), std(x)]);

N_point=101;
x_tran=min(x):(max(x)-min(x))/(N_point-1):max(x);

for ii=1:N_point
    y_tran(ii) = fitResult.a1 * exp(-((x_tran(ii) - fitResult.b1)/fitResult.c1)^2);
end

% % 绘制原始数据和拟合曲线
% plot(x, y, 'b.');
% hold on;
% plot(fitResult, 'r');
% hold off;
% 
% % 显示拟合结果
% disp(fitResult);
