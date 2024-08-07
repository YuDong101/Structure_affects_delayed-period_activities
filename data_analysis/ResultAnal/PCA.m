% 创建示例数据矩阵
X1 = F_Cc_avg_E;
X2 = F_Lens_avg_E;
X3 = F_EfromE;
X4 = F_EI_large;
X5 = F_All_Cycle;
Y = T_persistent_c';

% data = [X1, X2, X3, X4, X5, Y];
data = f_motifs;

% 计算数据的均值
meanData = mean(data);

% 对数据进行中心化
centeredData = data - meanData;

% 计算数据的协方差矩阵
covMatrix = cov(centeredData);

% 对协方差矩阵进行特征值分解
[eigenvectors, eigenvalues] = eig(covMatrix);

% 将特征值从大到小排序
[eigenvalues, idx] = sort(diag(eigenvalues), 'descend');
eigenvectors = eigenvectors(:, idx);

% 计算主成分得分
scores = centeredData * eigenvectors;

% 绘制主成分贡献率图
explainedVariance = eigenvalues / sum(eigenvalues);
figure;
bar(explainedVariance);
xlabel('主成分');
ylabel('贡献率');
title('主成分贡献率');

% 绘制主成分得分图
figure;
scatter(scores(:,1), scores(:,2));
xlabel('主成分1');
ylabel('主成分2');
title('主成分得分图');
%%
figure;
scatter(scores(:,10), Y);
xlabel('主成分1');
ylabel('Duration');
title('主成分得分图');


