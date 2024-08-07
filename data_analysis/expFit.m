close all
clc
clear X XX coefficients Y X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 R2 FStat pValue ErrorVariance initialGuess
at=find(T_persistent_c<=4000);

% X1 = normalize_data(P_Cc);
% X2 = normalize_data(P_lens);
% X3 = normalize_data(P_Dds_IE);
% X4 = normalize_data(Raito_richHub);
% X5 = normalize_data(P_Num_All_Cycle);

X1 = normalize_data(Cc_avg);
X2 = normalize_data(Cc_avg_E);
X3 = normalize_data(Lens_avg);
X4 = normalize_data(Lens_avg_E);
X5 = normalize_data(aa1);
X6 = normalize_data(aa2);
X7 = normalize_data(aa3);
X8 = normalize_data(aa4);
X9 = normalize_data(P_Dds_E_large);
X10 = normalize_data(P_Dds_I_large);
X11 = normalize_data(Num_All_Cycle_d);
X12 = normalize_data(Num_All_Cycle_E);
X13 = normalize_data(P_Cc);
X14 = normalize_data(P_lens);
X15 = normalize_data(P_Dds_IE);
X16 = normalize_data(Raito_richHub);
X17 = normalize_data(P_Num_All_Cycle);

y = T_persistent_c; % 先把这个转化成线性的

% 2. 定义多元指数回归模型
model = @(coeff, x) coeff(1) + exp(coeff(2)*x(:,1)) + exp(coeff(3)*x(:,2)) + ...
    exp(coeff(4)*x(:,3)) + exp(coeff(5)*x(:,4)) + exp(coeff(6)*x(:,5)) + ...
    exp(coeff(7)*x(:,6)) + exp(coeff(8)*x(:,7)) + exp(coeff(9)*x(:,8)) + ...
    exp(coeff(10)*x(:,9)) + exp(coeff(11)*x(:,10)) + exp(coeff(12)*x(:,11)) + ...
    exp(coeff(13)*x(:,12)) + exp(coeff(14)*x(:,13)) + exp(coeff(15)*x(:,14)) + ...
    exp(coeff(16)*x(:,15)) + exp(coeff(17)*x(:,16)) + exp(coeff(18)*x(:,17));

% model = @(coeff, x) coeff(1) + x(:,1).*coeff(2) + x(:,2).*coeff(3) + ...
%     x(:,3).*coeff(4) + x(:,4).*coeff(5) + x(:,5).*coeff(6) + ...
%     x(:,6).*coeff(7) + x(:,7).*coeff(8) + x(:,8).*coeff(9) + ...
%     x(:,9).*coeff(10) + x(:,10).*coeff(11) + x(:,11).*coeff(12) + ...
%     x(:,12).*coeff(13) + x(:,13).*coeff(14) + x(:,14).*coeff(15) + ...
%     x(:,15).*coeff(16) + x(:,16).*coeff(17) + x(:,17).*coeff(18);

% model = @(coeff, x) coeff(1) + exp(x(:,1).*coeff(2) + x(:,2).*coeff(3) + ...
%     x(:,3).*coeff(4) + x(:,4).*coeff(5) + x(:,5).*coeff(6) + ...
%     x(:,6).*coeff(7) + x(:,7).*coeff(8) + x(:,8).*coeff(9) + ...
%     x(:,9).*coeff(10) + x(:,10).*coeff(11) + x(:,11).*coeff(12) + ...
%     x(:,12).*coeff(13) + x(:,13).*coeff(14) + x(:,14).*coeff(15) + ...
%     x(:,15).*coeff(16) + x(:,16).*coeff(17) + x(:,17).*coeff(18));

% 4. 使用 lsqcurvefit 进行拟合
X = [X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15, X16, X17]; % 将五个自变量组合成一个矩阵
Y = y(at)'; % 调整 y 的格式
initialGuess(1:length(X(1,:))+1) = 0;
coefficients = lsqcurvefit(model, initialGuess, X(at,:), Y);
% coefficients = nlinfit(X(at,:), Y, model, initialGuess);

% 计算拟合后的因变量值
yFit = model(coefficients, X);

% 计算R² F p 残差
AAA=y';
SST = sum((AAA - mean(AAA)).^2); SSE = sum((AAA - yFit).^2); R2 = 1 - SSE / SST;
n = length(AAA); p = numel(coefficients) - 1; FStat = (SST - SSE) / p / (SSE / (n - p - 1));
pValue = 1 - fcdf(FStat, p, n-p-1);
ErrorVariance = SSE / (n-p-1);

% 显示统计信息
fprintf('R²（决定系数）: %.4f\n', R2);fprintf('F统计量: %.4f\n', FStat);fprintf('p值: %.4f\n', pValue);fprintf('估计的误差方差: %.4f\n', ErrorVariance);

[~, index_mlg] = sort(y);
%% semilogy
figure (301);
subplot(211), plot(1:length(index_mlg),y(index_mlg),1:length(index_mlg),yFit(index_mlg)); axis([1 length(at),0 6000]);  % semilogy
%%
function normalized_data = normalize_data(data)
    min_val = min(data);
    max_val = max(data);
    normalized_data = (data - min_val) / (max_val - min_val);
end