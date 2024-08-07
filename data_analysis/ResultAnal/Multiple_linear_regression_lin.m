
clc 
close all;
clear x y X Y X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 b_trail bb lamda_d best_lambda mse FitInfo B cv
rng(93)
at=find(T_persistent_c<=4000);
X1 = zscore(Cc_avg);
X2 = zscore(Cc_avg_E);
X3 = zscore(Lens_avg);
X4 = zscore(Lens_avg_E);
X5 = zscore(aaEE);
X6 = zscore(aaIE);
X7 = zscore(aaEI);
X8 = zscore(aaII);
X9 = zscore(P_Dds_E_large);
X10 = zscore(P_Dds_I_large);
X11 = zscore(Num_All_Cycle_d);
X12 = zscore(Num_All_Cycle_E);
X31 = zscore(P_Cc);
X32 = zscore(P_lens);
X33 = zscore(P_Dds_IE);
X34 = zscore(Raito_richHub);
X35 = zscore(P_Num_All_Cycle);
% X = ([X5, X6, X7, X8]); % , X1, X2, X31, X3, X4, X32, X5, X6, X7, X8, X33, X9, X10, X34, X11, X12, X35
X=[X1, X2, X31, X3, X4, X32, X5, X6, X7, X8, X33, X9, X10, X34, X11, X12, X35];

% X=[];
% P_motifs(isnan(P_motifs))=0;%% The number of motifs 13 is relatively small, and there may be a situation where f_motifs is 0. Exclude it
% for jj=1:13
%     X_Pmotifs(:,jj)=zscore(P_motifs(:,jj));
%     X_dmotifs(:,jj)=zscore(d_motifs(:,jj));
%     X_Emotifs(:,jj)=zscore(f_motifs_E(:,jj));
%     X_motifs(:,jj)=zscore(f_motifs(:,jj));
%     
%     X=[X X_Emotifs(at,jj)]; % X_Emotifs(:,jj) X_dmotifs(:,jj) X_Pmotifs(:,jj)
% end

Y=log(T_persistent_c');

x=[ones(size(X, 1), 1) X]; y=Y(at,:);

lambda_values = logspace(-3, 3, 100);
for ii=1:100
    
    mse_cv = zeros(size(lambda_values));
    for i = 1:length(lambda_values)
        lambda = lambda_values(i);
        mse_fold = zeros(5, 1);
        cv = cvpartition(length(at), 'KFold', 5);
        for j = 1:cv.NumTestSets
            train_idx = cv.training(j);
            test_idx = cv.test(j);
            X_train = x(train_idx, :);
            y_train = y(train_idx);
            X_test = x(test_idx, :);
            y_test = y(test_idx);
            [B(:,ii)] = ridge(y_train, X_train, lambda);
            y_pred = X_test * B(:,ii);
            mse_fold(j) = mean((y_test - y_pred).^2);
        end
        mse_cv(i) = mean(mse_fold);
    end

    [min_mse, min_idx] = min(mse_cv);
    best_lambda(ii) = lambda_values(min_idx);
    [b_trail(:,ii)] = ridge(y, x, best_lambda(ii));
    yPred = X * b_trail(2:end,ii) + b_trail(1,ii);
    mse(ii)=mean((y - yPred).^2);

    yMean = mean(Y);SS_total = sum((Y - yMean).^2);SS_residual = sum((Y - yPred).^2);
    R2(ii) = 1 - SS_residual / SS_total;
end

%%
[bb] = ridge(y, x, best_lambda(1));
ypred = X * bb(2:end) + bb(1);

[~, index_ranke] = sort(best_lambda); index_ranke=1;
figure(2);
subplot(411), semilogy(1:100,best_lambda(index_ranke),"ro");
subplot(412), plot(b_trail(2,index_ranke)',"ro");
subplot(413), semilogy(mse(index_ranke)',"ro");
meanb=mean(b_trail(2:end,:),2);
std_b=std(b_trail(2:end,:),0,2);
subplot(414),bar((meanb));
hold on;
errorbar(1:length(meanb), (meanb), std_b, 'k', 'linestyle', 'none');

figure(3);
mean_abs_bb=mean(abs(b_trail(2:end,:)),2);
std_abs_bb=std(abs(b_trail(2:end,:)),0,2);
subplot(311),bar((mean_abs_bb));
hold on;
errorbar(1:length(mean_abs_bb), (mean_abs_bb), std_abs_bb, 'k', 'linestyle', 'none');

meanb1=zscore(mean_abs_bb,1);
subplot(312),bar((meanb1));

log_b=log(abs(b_trail(2:end,:)));  log_b(isinf(log_b))=0;
meanb2=mean(log_b,2);std_b2=std(log_b,0,2);
subplot(313),bar(meanb2);
hold on;
errorbar(1:length(meanb2), meanb2, std_b2, 'k', 'linestyle', 'none');


[~, Rank_index] = sort(Y);
figure(4);
plot(exp(Y(Rank_index)), 'b', 'LineWidth', 2); axis([-inf inf,-inf inf]);  % semilogy
hold on;
plot(exp(ypred(Rank_index)), 'r--', 'LineWidth', 2); axis([-inf inf,-inf inf]);  % semilogy
grid on; hold off;


yMean = mean(Y);
SS_total = sum((Y - yMean).^2);
SS_residual = sum((Y - ypred).^2);
R2 = 1 - SS_residual / SS_total;


n = length(Y); 
p = size(x, 2); 
df_regression = p; 
df_residual = n - p - 1; 
F_statistic = (SS_total - SS_residual) / df_regression / (SS_residual / df_residual);

p_value = 1 - fcdf(F_statistic, df_regression, df_residual);

fprintf('R-squared (R^2): %.4f\n', R2);
fprintf('F statistic: %.4f\n', F_statistic);
fprintf('p value: %.4f\n', p_value);
